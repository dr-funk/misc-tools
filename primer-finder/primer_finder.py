#!/usr/bin/env python

import plotly.graph_objects as go
from Bio.Seq import Seq
import regex
import argparse
import time
from collections import defaultdict
import os
import sys

# parsing command line arguments
parser = argparse.ArgumentParser(description="a script to find primers matching your sequence of interest")
parser.add_argument("-s", "--sequence", help="the sequence you want to analyze. provide a path to a fasta file or just paste the sequence here", type=str, metavar='<path/ATCG>')
parser.add_argument("-o", "--outdir", help="folder in which the results should be saved. if not provided results will be displayed in browser", type=str, metavar='<path>')
parser.add_argument("-p", "--primers", help='path to the text file containing primers to be used', metavar='<path>')
parser.add_argument('-e','--edge', help="no mismatches allowed within N nucleotides of the 3' end", type=int, default=5, metavar='N')
parser.add_argument("-mm", "--mismatches", help='how many mismatches are allowed', default=1, type=int, metavar='N')
parser.add_argument('-d','--degenerate', help='allow at most this many degnerate positions', default=3, type=int, metavar='N')
parser.add_argument('-to','--textonly', help='generate only the text output file', action='store_true')
parser.add_argument('-go', '--graphonly', help='generate only the interactive graph output', action='store_true')

args=parser.parse_args()

# some primers contain chars we don't want
# remove all spaces, replace U with T
# sometime there's illegal chars, replace with N if not in IUPAC
def fix_primer(primer):
    allowed_chars='AGCTRYMKSWHBVDN'
    primer=primer.upper()
    primer.replace('T','U')
    primer.replace(' ','')
    primer=''.join((c if c in allowed_chars else 'N' for c in primer))
    return primer

def deambiguoize(primer):
    for c in primer:
        if c not in ('A','T','C','G'):
            # at least one ambiguous char
            primer=''.join((c if c not in degenerate.keys() else degenerate[c] for c in primer))
            return primer
    else:
        # no ambiguous
        return primer
    
# I want to skip taqman probes, these should be recognizable by having dyes listed in field 15
dyes=("5 FAM 3'BHQ1", 'Cy5-BHQ3', 'FAM-BHQ', 'VIC-QSY', 'Atto425-DDQ', 'DRAGONFLY-BHQ2', 'FAM', 'FAM - BHQ1', 'FAM-MGB', '6-FAM-MGB-Eclipse', 'JUN-QSY', 'NED', 'VIC-TAMRA', 'VIC', 'Cy5', 'ABY-QSY', '6-FAM', "5' FAM 3'BHQ-1", 'FAM-TAMRA', 'EE-1775a (PCR room)', 'DFO-BHQ2', 'JOE', 'Biotine', 'ATTO425-BHQ1', 'FAM-NFQ', 'DFO', 'PET', 'Ee-1775a (PCR room)', 'Atto425')

# read in primers, doing it like this since some lines are incomplete
if not args.primers:
    print('No primer file path provided, attempting to locate in script directory')
    try:
        primers=defaultdict(list)
        with open(os.path.dirname(os.path.realpath(__file__))+os.sep+'RF-primers.txt','r') as f:
            for line in f:
                l=line.strip().split('\t')
                if len(l)>=16 and l[15] not in dyes:
                    primers[fix_primer(l[7].strip('"'))].append(l[1].strip('"'))
        print(f'loaded primers from {sys.argv[0].rsplit(os.sep,1)[0]}{os.sep}RF-primers.txt')
    except:
        sys.exit('Error: no primer file provided! please use -p to tell me where the primer data is.')
else:
    primers=defaultdict(list)
    with open(args.primers,'r') as f:
        for line in f:
            l=line.strip().split('\t')
            if len(l)>=16 and l[15] not in dyes:
                primers[fix_primer(l[7].strip('"'))].append(l[1].strip('"'))


seq_input = args.sequence
# check if this a path to a fasta file
if not args.sequence:
    seq_id=''
    target_seq=input('paste sequence now:')
else:
    if os.path.isfile(args.sequence):
        with open(args.sequence, 'r') as f:
            seq_id=f.readline().strip()[1:].replace(os.sep, '-')+'_'
            target_seq=''.join(tuple(line.strip() for line in f)) # concatenate all remaining lines
            target_seq=target_seq.upper().replace('U','T')
    else:
        seq_id=''
        target_seq=args.sequence.upper().replace('U','T')

dist_thrshld = args.mismatches
# build the color dictionary by interpolation of black to red
colors={m:f'#{int((m/args.mismatches)*255):02x}0000' for m in range(args.mismatches+1)} if args.mismatches != 0 else {0:'#000000'}
end_thrshld = args.edge
degen_thrshld = args.degenerate

degenerate={
    'R':'[AG]',
    'Y':'[CT]',
    'M':'[AC]',
    'K':'[GT]',
    'S':'[CG]',
    'W':'[AT]',
    'H':'[ACT]',
    'B':'[CGT]',
    'V':'[ACG]',
    'D':'[AGT]',
    'N':'[ACGT]',
}

t0 = time.time()

primer_matches=[]
# iterate through primers and align
for n, (primer_seq, primer_nbs) in enumerate(primers.items()):
    primer_seq_rc=str(Seq(primer_seq).reverse_complement())

    if (n+1) % 1000 == 0:
        print(f'\raligning primer {n+1}/{len(primers)}', end='')

    # check how many degenerate characters we have
    if sum(1 for c in primer_seq if c not in ('ATCG')) > degen_thrshld:
        continue

    # fuzzy regex matching forward and rev
    if dist_thrshld!=0:
        # building patterns: 2 parts: 
        #   the 3' part has to match perfectly
        #   the 5' part can have mismatches
        pattern_fwd = regex.compile('(?b)'+ # look for best match
                                    f'(?:{deambiguoize(primer_seq[:-end_thrshld])}){{s<={dist_thrshld}}}'+ # fuzzy match 5' 
                                    f'(?:{deambiguoize(primer_seq[-end_thrshld:])})' # perfect match 3'
                                    )
        pattern_rev = regex.compile('(?b)'+ # look for best match
                                    f'(?:{deambiguoize(primer_seq_rc[:-end_thrshld])}){{s<={dist_thrshld}}}'+ # fuzzy match 5' 
                                    f'(?:{deambiguoize(primer_seq_rc[-end_thrshld:])})' # perfect match 3'
                                    )
    else:
        pattern_fwd = regex.compile(deambiguoize(primer_seq))
        pattern_rev = regex.compile(deambiguoize(primer_seq_rc))

    # find all matches of the primer, with details in match object
    fwd_matches = tuple(pattern_fwd.finditer(target_seq, overlapped=True))
    rev_matches = tuple(pattern_rev.finditer(target_seq, overlapped=True))

    if dist_thrshld>0 and len(fwd_matches+rev_matches)>1:
        # filter the matches to only keep matches with the lowest number of errors
        # else primers that match once perfectly and once with errors will be discarded
        min_mm=min(len(match.fuzzy_changes[0]) for match in fwd_matches+rev_matches)

        fwd_matches=tuple(match for match in fwd_matches if len(match.fuzzy_changes[0]) == min_mm)
        rev_matches=tuple(match for match in rev_matches if len(match.fuzzy_changes[0]) == min_mm)

    if len(fwd_matches+rev_matches) == 1:
        # there is only one single valid match across orientations
        if not rev_matches:
            ori, match = 'fwd', fwd_matches[0]
            changes=','.join((target_seq[i]+str(i)+primer_seq[i-match.span(0)[0]] for i in match.fuzzy_changes[0]))
        else:
            ori, match = 'rev', rev_matches[0]
            changes=','.join((target_seq[i]+str(i)+primer_seq_rc[i-match.span(0)[0]] for i in match.fuzzy_changes[0]))
    else:
        # either no or too many matches
        continue

    # if we're here, we found exaclty one valid match, add it to the list

    primer_matches.append((','.join(primer_nbs), 
                           ori,
                           match.span(0)[0]+1,
                           match.span(0)[1]+1,
                           changes if changes else 'none',
                           primer_seq if ori=='fwd' else primer_seq_rc,
                           match.group(0)
                            ))
print(f'\raligning primer {n+1}/{len(primers)}', end='')
print(f'\ndone, took {(time.time()-t0):.2f} seconds, found {len(primer_matches)} hits')

# sort by match start
primer_matches.sort(key=lambda x: x[2])

if not args.graphonly and args.outdir:
    # write to file
    with open(os.path.join(args.outdir,f'{seq_id}out.tsv'), 'w') as fout:
        for match in primer_matches:
            fout.write('\t'.join(map(str, match)))
            fout.write('\n')

# make graph
layers_pos=defaultdict(list)
layers_pos[0]=[0]
layers_marker=defaultdict(list)
layers_hover=defaultdict(list)
layers_color=defaultdict(list)
ori_markers={'fwd':('line-ns','arrow-bar-right'),
             'rev':('arrow-bar-left','line-ns')}

# go through all the primer hits, collect data
for match in primer_matches:
    primer_nbs, ori, start, end, muts, primer_seq, match_seq = match

    # build a formatted alignment
    if muts=='none':
        mut_pos=tuple()
        aln='|'*(end-start)
    else:
        aln=''
        mut_pos = tuple(int(mut[1:-1])+1 for mut in muts.split(','))
        for i in range(start, end):
            if i in mut_pos:
                aln+=' '
            else:
                aln+='|'

    # fit the primer into a free part of a layer
    for layer in layers_pos.keys():
        if layers_pos[layer][-1]<start:
            free_layer=layer
            break
    else:
        free_layer=layer+0.5

    layers_pos[free_layer].extend((start, end))
    layers_marker[free_layer].extend((ori_markers[ori][0],ori_markers[ori][1]))
    layers_hover[free_layer].extend(((primer_nbs,'forward' if ori == 'fwd' else 'reverse',start,end,muts, primer_seq, match_seq,aln),
                                     (primer_nbs,'forward' if ori == 'fwd' else 'reverse',start,end,muts, primer_seq, match_seq,aln)))
    layers_color[free_layer].extend((colors[len(mut_pos)], colors[len(mut_pos)]))
   
if not args.textonly:
    fig=go.Figure()
    layers_pos[0]=layers_pos[0][1:]
    for layer in layers_pos.keys():

        for x0,x1,color in zip(layers_pos[layer][0::2], layers_pos[layer][1::2],layers_color[layer][0::2]):
            fig.add_shape(type='line',
                    x0=x0,
                    y0=layer,
                    x1=x1,
                    y1=layer,
                    line=dict(color=color,),
                    xref='x',
                    yref='y',
                    line_width=2,
                    opacity=1,
                    )
        fig.add_trace(go.Scatter(
                            x=layers_pos[layer],
                            y=[layer]*len(layers_pos[layer]),
                            mode='markers',
                            customdata=layers_hover[layer],
                            hovertemplate='<b>%{customdata[0]}</b><br>in %{customdata[1]} orientation %{customdata[2]}-%{customdata[3]}<br>mutations: %{customdata[4]}<br>'+
                                            '%{customdata[6]}<br>%{customdata[7]}<br>%{customdata[5]}</span><extra></extra>',
                            marker=dict(
                                    symbol=layers_marker[layer],
                                    size=15,
                                    color=layers_color[layer],
                                    line=dict(width=2, color=layers_color[layer])
                                    ),

                                )
                            )

    fig.update_layout(hoverlabel=dict(bgcolor='white',
                                        font_family='Courier'),
                    showlegend=False,
                    template='simple_white',
                    )
    fig.update_xaxes(range=[-1, len(target_seq)+2])
    fig.update_yaxes(visible=False)

    if args.outdir:
        fig.write_html(os.path.join(args.outdir, f'{seq_id}out.html'), include_plotlyjs='cdn')
    else:
        fig.show()

