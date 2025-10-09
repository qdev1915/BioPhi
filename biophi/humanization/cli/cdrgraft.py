from functools import partial
from multiprocessing import Pool

import click
from Bio import SeqIO
import pandas as pd
from biophi.common.utils.formatting import logo
from biophi.common.utils.io import parse_antibody_files, write_sheets
from biophi.common.utils.seq import iterate_fasta
from biophi.humanization.cli.oasis import show_unpaired_warning
from biophi.humanization.methods.humanness import OASisParams, get_antibody_humanness
from biophi.humanization.methods.humanization import (
    humanize_antibody, 
    CDRGraftingHumanizationParams, 
    HumanizationParams
)
from abnumber import Chain, ChainParseError, SUPPORTED_CDR_DEFINITIONS, SUPPORTED_SCHEMES
import os
import sys
from tqdm import tqdm


@click.command()
@click.argument('inputs', required=False, nargs=-1)
@click.option('--output', required=False, help='Output directory path. With --fasta-only, output FASTA file path.')
@click.option('--fasta-only', is_flag=True, default=False, type=bool, help='Output only a FASTA file with humanized sequences (speeds up processing)')
@click.option('--oasis-db', required=False, help='OAS peptide database connection string (required to run OASis)')
@click.option('--scheme', default=HumanizationParams.scheme, help=f'Numbering scheme: one of {", ".join(SUPPORTED_SCHEMES)}')
@click.option('--cdr-definition', default=HumanizationParams.cdr_definition, help=f'CDR definition: one of {", ".join(SUPPORTED_CDR_DEFINITIONS)}')
@click.option('--heavy-v-germline', default='auto', help='Heavy chain V germline gene (auto for automatic selection)')
@click.option('--light-v-germline', default='auto', help='Light chain V germline gene (auto for automatic selection)')
@click.option('--backmutate-vernier/--no-backmutate-vernier', default=False, type=bool, help='Backmutate Vernier zone residues to parental (default: disabled)')
@click.option('--sapiens-iterations', type=int, default=0, help='Additional Sapiens iterations after CDR grafting')
@click.option('--limit', required=False, metavar='N', type=int, help='Process only first N records')
def cdrgraft(inputs, output, fasta_only, scheme, cdr_definition, heavy_v_germline, light_v_germline, 
             backmutate_vernier, sapiens_iterations, limit, oasis_db):
    """CDR Grafting: Humanize antibodies by grafting CDRs onto human germline frameworks.

    CDR Grafting transplants the CDR loops from the input antibody onto the most similar
    human germline framework, optionally with Vernier zone backmutations and Sapiens refinement.

    EXAMPLES:

        \\b
        # Graft CDRs onto auto-selected human germlines, print to standard output
        biophi cdrgraft input.fa

        \\b
        # CDR graft with Vernier zone backmutations (default)
        biophi cdrgraft input.fa --fasta-only --output humanized.fa

        \\b
        # CDR graft followed by 2 iterations of Sapiens refinement
        biophi cdrgraft input.fa --sapiens-iterations 2 --output ./report/ \\\\
          --oasis-db /path/to/OASis_9mers_v1.db

        \\b
        # Specify custom germline genes
        biophi cdrgraft input.fa --heavy-v-germline IGHV3-23*01 \\\\
          --light-v-germline IGKV1-39*01 --output ./report/

    INPUTS: Input FASTA file path(s). If not provided, creates an interactive session.
    """

    logo('''    ____ ____  ____     ____            __ _   
   / ___|  _ \\|  _ \\   / ___|_ __ __ _/ _| |_ 
  | |   | | | | |_) | | |  _| '__/ _` | |_| __|
  | |___| |_| |  _ <  | |_| | | | (_| |  _| |_ 
   \\____|____/|_| \\_\\  \\____|_|  \\__,_|_|  \\__|
                                                  ''')

    click.echo(f'Settings:', err=True)
    click.echo(f'- Humanization method: CDR Grafting', err=True)
    click.echo(f'- Numbering scheme: {scheme}', err=True)
    click.echo(f'- CDR definition: {cdr_definition}', err=True)
    click.echo(f'- Heavy V germline: {heavy_v_germline}', err=True)
    click.echo(f'- Light V germline: {light_v_germline}', err=True)
    click.echo(f'- Backmutate Vernier zone: {"Yes" if backmutate_vernier else "No"}', err=True)
    if sapiens_iterations > 0:
        click.echo(f'- Sapiens refinement iterations: {sapiens_iterations}', err=True)
    click.echo(f'', err=True)

    humanization_params = CDRGraftingHumanizationParams(
        scheme=scheme,
        cdr_definition=cdr_definition,
        heavy_v_germline=heavy_v_germline,
        light_v_germline=light_v_germline,
        backmutate_vernier=backmutate_vernier,
        sapiens_iterations=sapiens_iterations
    )
    
    oasis_params = OASisParams(
        oasis_db_path=oasis_db,
        min_fraction_subjects=0.10
    ) if oasis_db else None

    if inputs:
        if fasta_only:
            return cdrgraft_fasta_only(
                inputs,
                output,
                humanization_params=humanization_params,
                limit=limit
            )
        else:
            return cdrgraft_full_report(
                inputs,
                output,
                humanization_params=humanization_params,
                oasis_params=oasis_params,
                limit=limit
            )
    else:
        # Interactive mode
        return cdrgraft_interactive(
            humanization_params=humanization_params,
            oasis_params=oasis_params
        )


def cdrgraft_fasta_only(inputs, output, humanization_params, limit=None):
    """Process FASTA files and output only humanized sequences."""
    click.echo('Reading input files...', err=True)
    records = list(iterate_fasta(inputs))
    
    if limit:
        records = records[:limit]
    
    click.echo(f'Processing {len(records)} sequences...', err=True)
    
    humanized_records = []
    
    for record in tqdm(records, desc='Humanizing', file=sys.stderr):
        try:
            chain = Chain(str(record.seq), scheme=humanization_params.scheme, cdr_definition=humanization_params.cdr_definition)
            chain.name = record.id
            
            # Determine if heavy or light chain
            if chain.is_heavy_chain():
                result = humanize_antibody(vh=chain, vl=None, params=humanization_params)
                humanized_chain = result.vh.humanized_chain if result.vh else None
            else:
                result = humanize_antibody(vh=None, vl=chain, params=humanization_params)
                humanized_chain = result.vl.humanized_chain if result.vl else None
            
            if humanized_chain:
                # Create description
                chain_type = 'VH' if humanized_chain.is_heavy_chain() else 'VL'
                method_desc = f'CDR_Grafted_{humanization_params.cdr_definition}_'
                if humanization_params.backmutate_vernier:
                    method_desc += 'Vernier_'
                if humanization_params.sapiens_iterations > 0:
                    method_desc += f'Sapiens_{humanization_params.sapiens_iterations}iter_'
                
                germline = humanized_chain.v_gene if hasattr(humanized_chain, 'v_gene') else 'auto'
                description = f'{record.id} {chain_type} (Humanized {record.id} {method_desc}{germline} BioPhi)'
                
                from Bio.SeqRecord import SeqRecord
                humanized_record = SeqRecord(
                    seq=humanized_chain.seq,
                    id=record.id,
                    description=description
                )
                humanized_records.append(humanized_record)
        
        except ChainParseError as e:
            click.echo(f'Warning: Could not parse {record.id}: {e}', err=True)
            continue
        except Exception as e:
            click.echo(f'Error processing {record.id}: {e}', err=True)
            continue
    
    # Output results
    if output:
        click.echo(f'Writing output to {output}...', err=True)
        with open(output, 'w') as f:
            SeqIO.write(humanized_records, f, 'fasta')
    else:
        # Print to stdout
        SeqIO.write(humanized_records, sys.stdout, 'fasta')
    
    click.echo(f'Completed: {len(humanized_records)} sequences humanized', err=True)


def cdrgraft_full_report(inputs, output, humanization_params, oasis_params=None, limit=None):
    """Process FASTA files and generate full report with OASis analysis."""
    click.echo('Reading input files...', err=True)
    
    # Read FASTA records
    records = list(iterate_fasta(inputs))
    
    if limit:
        records = records[:limit]
    
    if not records:
        click.echo('No valid sequences found!', err=True)
        return
    
    click.echo(f'Processing {len(records)} sequences...', err=True)
    
    results = []
    humanized_records = []
    
    for record in tqdm(records, desc='Humanizing', file=sys.stderr):
        try:
            chain = Chain(str(record.seq), scheme=humanization_params.scheme, cdr_definition=humanization_params.cdr_definition)
            chain.name = record.id
            
            # Determine if heavy or light chain
            if chain.is_heavy_chain():
                vh_chain = chain
                vl_chain = None
                result = humanize_antibody(vh=vh_chain, vl=None, params=humanization_params)
            else:
                vh_chain = None
                vl_chain = chain
                result = humanize_antibody(vh=None, vl=vl_chain, params=humanization_params)
            
            # Get humanness scores if OASis DB is available
            parental_humanness = None
            humanized_humanness = None
            
            if oasis_params:
                try:
                    parental_humanness = get_antibody_humanness(
                        vh=vh_chain,
                        vl=vl_chain,
                        params=oasis_params
                    )
                    humanized_humanness = get_antibody_humanness(
                        vh=result.vh.humanized_chain if result.vh else None,
                        vl=result.vl.humanized_chain if result.vl else None,
                        params=oasis_params
                    )
                except Exception as e:
                    click.echo(f'Warning: Could not compute OASis scores for {record.id}: {e}', err=True)
            
            # Store result
            results.append({
                'name': record.id,
                'humanization': result,
                'parental_humanness': parental_humanness,
                'humanized_humanness': humanized_humanness
            })
            
            # Create FASTA records
            for chain_result in [result.vh, result.vl]:
                if chain_result:
                    humanized_chain = chain_result.humanized_chain
                    chain_type = 'VH' if humanized_chain.is_heavy_chain() else 'VL'
                    
                    method_desc = humanization_params.get_export_name()
                    description = f'{record.id} {chain_type} (Humanized {record.id} {method_desc}BioPhi)'
                    
                    from Bio.SeqRecord import SeqRecord
                    humanized_record = SeqRecord(
                        seq=humanized_chain.seq,
                        id=f'{record.id}_{chain_type}',
                        description=description
                    )
                    humanized_records.append(humanized_record)
        
        except Exception as e:
            click.echo(f'Error processing {record.id}: {e}', err=True)
            continue
    
    # Output results
    if not output:
        # Print to stdout
        SeqIO.write(humanized_records, sys.stdout, 'fasta')
    else:
        # Create output directory
        os.makedirs(output, exist_ok=True)
        
        # Write humanized FASTA
        fasta_path = os.path.join(output, 'humanized.fa')
        click.echo(f'Writing humanized sequences to {fasta_path}...', err=True)
        with open(fasta_path, 'w') as f:
            SeqIO.write(humanized_records, f, 'fasta')
        
        # Write alignments file
        if results:
            alignments_path = os.path.join(output, 'alignments.txt')
            click.echo(f'Writing alignments to {alignments_path}...', err=True)
            alignments = [res['humanization'].get_alignment_string() for res in results]
            with open(alignments_path, 'w') as f:
                f.write('\n\n'.join(alignments))
        
        # Write Excel report if OASis was run
        if oasis_params and results:
            xlsx_path = os.path.join(output, 'CDRGraft_report.xlsx')
            click.echo(f'Writing report to {xlsx_path}...', err=True)
            
            sheets = {}
            
            # Overview sheet
            overview_data = []
            for res in results:
                row = {'Antibody': res['name']}
                
                if res['humanized_humanness']:
                    hum = res['humanized_humanness']
                    min_frac = oasis_params.min_fraction_subjects if oasis_params else 0.01
                    row['OASis_Identity'] = hum.get_oasis_identity(min_frac)
                    row['OASis_Percentile'] = hum.get_oasis_percentile(min_frac)
                    if hum.vh:
                        row['Heavy_OASis_Identity'] = hum.vh.get_oasis_identity(min_frac)
                        # row['Heavy_Germline'] = hum.vh.v_gene  # Not available in ChainHumanness
                    if hum.vl:
                        row['Light_OASis_Identity'] = hum.vl.get_oasis_identity(min_frac)
                        # row['Light_Germline'] = hum.vl.v_gene  # Not available in ChainHumanness
                
                num_mutations = 0
                if res['humanization'].vh:
                    num_mutations += res['humanization'].vh.num_mutations()
                if res['humanization'].vl:
                    num_mutations += res['humanization'].vl.num_mutations()
                row['Num_Mutations'] = num_mutations
                
                overview_data.append(row)
            
            sheets['Overview'] = pd.DataFrame(overview_data)
            
            # Write to Excel
            write_sheets(sheets, xlsx_path)
        
        click.echo(f'Completed! Output saved to {output}', err=True)
    
    click.echo(f'Successfully humanized {len(results)} {"antibody" if len(results) == 1 else "antibodies"}', err=True)


def cdrgraft_interactive(humanization_params, oasis_params=None):
    """Interactive mode for single antibody humanization."""
    click.echo('Interactive mode - Enter antibody sequences:', err=True)
    click.echo('(Paste VH and/or VL sequences, press Enter twice when done)', err=True)
    
    sequences = []
    while True:
        try:
            line = input()
            if not line.strip():
                if sequences:
                    break
                continue
            sequences.append(line.strip())
        except EOFError:
            break
    
    if not sequences:
        click.echo('No sequences provided', err=True)
        return
    
    # Parse sequences
    vh_chain = None
    vl_chain = None
    
    for seq in sequences:
        if seq.startswith('>'):
            continue
        try:
            chain = Chain(seq, scheme=humanization_params.scheme, cdr_definition=humanization_params.cdr_definition)
            if chain.is_heavy_chain():
                vh_chain = chain
                vh_chain.name = 'VH'
            else:
                vl_chain = chain
                vl_chain.name = 'VL'
        except ChainParseError:
            continue
    
    if not vh_chain and not vl_chain:
        click.echo('Could not parse any valid antibody sequences', err=True)
        return
    
    # Humanize
    click.echo('\nHumanizing...', err=True)
    result = humanize_antibody(vh=vh_chain, vl=vl_chain, params=humanization_params)
    
    # Print alignment
    click.echo('\n' + result.get_alignment_string())
    
    # Print humanness if available
    if oasis_params:
        try:
            humanness = get_antibody_humanness(
                vh=result.vh.humanized_chain if result.vh else None,
                vl=result.vl.humanized_chain if result.vl else None,
                params=oasis_params
            )
            click.echo(f'\nOASis Identity: {humanness.get_oasis_identity(oasis_params.min_fraction_subjects):.2%}')
            click.echo(f'OASis Percentile: {humanness.get_oasis_percentile(oasis_params.min_fraction_subjects):.2%}')
        except Exception as e:
            click.echo(f'\nCould not compute OASis scores: {e}', err=True)


if __name__ == '__main__':
    cdrgraft()
