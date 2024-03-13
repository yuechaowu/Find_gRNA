#! /usr/bin/env python

import sys
from os import path,makedirs,cpu_count
import gc
import regex as re
from pyfastx import Fasta
import datetime
import subprocess
from multiprocessing import Pool
import argparse


def reverse_complementary(seq, reverse=True, compelementary=True):
    """

    :param seq:
    :param reverse:

        False:

            original:  ACAACCCACTAAACCATTACTGG
            converted: ACAACCCACTAAACCATTACTGG

        True:

            original:  ACAACCCACTAAACCATTACTGG
            converted: GGTCATTACCAAATCACCCAACA
    :param compelementary:

        False:

            original:  ACAACCCACTAAACCATTACTGG
            converted: ACAACCCACTAAACCATTACTGG

        True:

            original:  ACAACCCACTAAACCATTACTGG
            converted: TGTTGGGTGATTTGGTAATGACC

    :return:
    """

    str_trans = str.maketrans("ATGCRYSWKMBDHV", "TACGYRWSMKVHDB")

    if reverse and not compelementary:                                   # Reverse

        rev_seq = str(seq)[::-1]

    if compelementary and not reverse:                                   # compelementary

        rev_seq = str(seq).translate(str_trans)

    if reverse and compelementary:                                       # Only do reverse

        rev_seq = str(seq).translate(str_trans)[::-1]

    if not reverse and not compelementary:

        rev_seq = str(seq)

    return rev_seq


def calc_homopolymers(sequence, l=5):

    if 'A'*l in sequence or 'T'*l in sequence or 'C'*l in sequence or 'G'*l in sequence:

        homopolymers = True

    else:

        homopolymers = False

    return homopolymers


def calc_GC(sequence):

    G_count = sequence.upper().count('G')

    C_count = sequence.upper().count('C')

    GC_content = float((G_count + C_count) / len(sequence))

    return GC_content


def find_putative_seq(seqfile):
    """

    :param seqfile:
        sequence.fasta

    :param length:
        NGG:  20
        TTTV: 23

    :param pam:
        SpCas9 from Streptococcus pyogenes: 5'-NGG-3'
        SpCas9 from Streptococcus pyogenes: 5'-NGG-3' for target, 5'-NRG-3' for off-target (R = A or G)
        StCas9 from Streptococcus thermophilus: 5'-NNAGAAW-3' (W = A or T)
        NmCas9 from Neisseria meningitidis: 5'-NNNNGMTT-3' (M = A or C)
        SaCas9 from Staphylococcus aureus: 5'-NNGRRT-'3 (R=A or G)
        CjCas9 from Campylobacter jejuni: 5'-NNNVRYAC-3' (V = G or C or A, R = A or G, Y = C or T)
        CjCas9 from Campylobacter jejuni: 5'-NNNNRYAC-3' (R = A or G, Y = C or T)
        AsCpf1 from Acidaminococcus or LbCpf1 from Lachnospiraceae: 5'-TTTN-3'
        AsCpf1 from Acidaminococcus or LbCpf1 from Lachnospiraceae: 5'-TTTV-3' (V = G or C or A)
        SpCas9 from Streptococcus pasteurianus: 5'-NNGTGA-3'
        FnCpf1 from Francisella: 5'-TTN-3'
        SaCas9 from Staphylococcus aureus: 5'-NNNRRT-'3 (R=A or G)
        FnCpf1 from Francisella: 5'-KYTV-3'
        VRER SpCas9 from Streptococcus pyogenes: 5'-NGCG-3'
        VQR SpCas9 from Streptococcus pyogenes: 5'-NGA-3'
        XCas9 3.7 (TLIKDIV SpCas9) from Streptococcus pyogenes: 5'-NGT-3'
        XCas9 3.7 (TLIKDIV SpCas9) from Streptococcus pyogenes: 5'-NG-3'

    :param isreversed:

        False:
            NNNNNNNNNNNNNNNNNNNNNGG

        True:
            TTTVNNNNNNNNNNNNNNNNNNNNNNN

    :return:
    """

    fasta = Fasta(seqfile)

    if isreversed:                                       # PAM followed by gRNA (PAM + gRNA)

        p = pam + length*'N'

        repl = "^\w{" + str(len(pam)) + '}'

    else:                                                # gRNA followed by PAM (gRNA + PAM)

        p = length*'N' + pam

        repl = "\w{" + str(len(pam)) + '}$'

    p_rev = reverse_complementary(p, True, False)        # reverse gRNA pattern

    pattern = p.replace('N', '[AGTC]').replace('R', '[AG]').replace('W', '[AT]').replace('M', '[AC]').replace('Y', '[CT]').replace('S', '[GC]').replace('K', '[GT]').replace('B', '[CGT]').replace('D', '[AGT]').replace('H', '[ACT]').replace('V', '[ACG]')

    pattern_rev = p_rev.replace('N', '[AGTC]').replace('R', '[AG]').replace('W', '[AT]').replace('M', '[AC]').replace('Y', '[CT]').replace('S', '[GC]').replace('K', '[GT]').replace('B', '[CGT]').replace('D', '[AGT]').replace('H', '[ACT]').replace('V', '[ACG]')

    match = re.compile(pattern)

    match_rev = re.compile(pattern_rev)

    sequences = set()

    ot_seq = set()

    outio = open(work_dir+'/gRNA_candidates.txt', 'w')

    for target_seq in fasta:

        gene = target_seq.name

        # target_seq = fasta[gene][:]
        target_seq = target_seq.seq

        for seq in match.finditer(target_seq, overlapped=True):           # Find 5' -> 3' gRNA

            chromosome = gene

            sequence = seq.group()

            sgrna_no_pam = re.sub(repl, '', sequence)

            remove = calc_homopolymers(sgrna_no_pam, 5)

            if not remove:

                ot_seq.add(re.sub(repl, 'N'*len(pam), sequence))

            if isreversed:

                start = str(seq.span()[0]+len(pam)+1)

                end = str(seq.span()[1])

            else:

                start = str(seq.span()[0]+1)

                end = str(seq.span()[1]-len(pam))

            span = chromosome + ':' + str(seq.span()[0]+1) + '-' + str(seq.span()[1])

            sequences.add(chromosome + '\t' + start + '\t' + end + '\t' + sgrna_no_pam + '\t' + span + '\t' + sequence + '\t+')

        for seq in match_rev.finditer(reverse_complementary(target_seq, False, True), overlapped=True):       # Find 3' -> 5' gRNA

            chromosome = gene

            sequence = seq.group()

            sgrna_no_pam = re.sub(repl, '', reverse_complementary(sequence, True, False))

            remove = calc_homopolymers(sgrna_no_pam, 5)

            if not remove:

                ot_seq.add(re.sub(repl, 'N'*len(pam), reverse_complementary(sequence, True, False)))

            if isreversed:

                start = str(seq.span()[0]+1)

                end = str(seq.span()[1]-len(pam))

            else:

                start = str(seq.span()[0]+len(pam)+1)

                end = str(seq.span()[1])

            span = chromosome + ':' + str(seq.span()[0]+1) + '-' + str(seq.span()[1])

            sequences.add(chromosome + '\t' + start + '\t' + end + '\t' + sgrna_no_pam + '\t' + span + '\t' + sequence + '\t-')

    for seq in sequences:

        print(seq, file = outio)

    outio.close()

    off_out = open(work_dir+'/offinder_configure.txt', 'w')

    print(path.abspath(path.expanduser(offdir)), file=off_out)
    print(p, file=off_out)

    for seqinfo in ot_seq:

        print(seqinfo, mismatch, file=off_out)

    off_out.close()


def cas_offinder(offinder, device):

    check_types = subprocess.Popen(offinder, shell=True, stdout=subprocess.PIPE).stdout.read().decode()

    if 'Type: GPU' in check_types:

        dev = 'G'

    elif 'Type: CPU' in check_types:

        dev = 'C'

    else:

        print('No device support cas-offinder, please check if the OpenCL/OpenGL or GPU drivers installed.')

        return 0

    if device:

        dev = device

    print('Use %sPU to run cas-offinder program.' % dev)

    command = offinder + ' ' + work_dir+'/offinder_configure.txt ' + dev + ' '+work_dir+'/off-target_sites.txt'

    print('cas-offinder command: ', command, '', sep='\n')

    subprocess.call(command, shell=True)

    off_targets = dict()

    with open(work_dir+'/off-target_sites.txt', 'r') as offinfo:

        for info in offinfo:

            sgRNA_no_pam = info.rstrip().split('\t')[0].replace('N'*len(pam), '')
            mismatch_info = info.rstrip().split('\t')[-1]

            if sgRNA_no_pam not in off_targets:

                off_targets[sgRNA_no_pam] = dict()

            if mismatch_info not in off_targets[sgRNA_no_pam]:

                off_targets[sgRNA_no_pam][mismatch_info] = 1

            else:

                off_targets[sgRNA_no_pam][mismatch_info] += 1

    offinfo.close()

    outf = open(outfile, 'w')

    with open(work_dir+'/gRNA_candidates.txt', 'r') as all_info:

        for line in all_info:

            if line.startswith('#'):

                continue

            info = line.rstrip().split('\t')
            grna_no_pam = info[3]
            grna_strand = info[6]

            GC_content = calc_GC(grna_no_pam)

            name = dict()

            for i in range(mismatch+1):

                miscount = 'm'+str(i)

                if str(i) in off_targets[grna_no_pam]:

                    name[miscount] = off_targets[grna_no_pam][str(i)]
                else:

                    name[miscount] = 0
            
            print('\t'.join(info), GC_content, sep='\t', end='', file=outf)

            for miscount in range(mismatch, -1, -1):

                print("\t"+str(name['m'+str(miscount)]), end='', file=outf)

            print(file=outf)
    
    outf.close()

    all_info.close()

    return 1


def get_options():

    parser = argparse.ArgumentParser(description="Program for finding position of CRISPR gRNA in specific sequences", prog='CRISPRPosFinder', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="Example:\n"
                                            "  python crisprposfinder.py -s genome.fa -l 20 -p NGG -m 3 --offind \\\n"
                                            "                            --offpath /path_to_program/cas-offinder \\\n"
                                            "                            -o results/cas_output.txt")

    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+'v1.0.1')

    parser.add_argument('-s', dest='seqfile', required=True, type=str, help='Input sequence fasta file')

    parser.add_argument('-l', dest='length', required=True, type=int, help='Length of guide RNA without PAM')

    parser.add_argument('-p', dest='pam', required=True, type=str, help='CRISPR PAM sequence')

    parser.add_argument('-m', dest='mismatch', type=int, default=3, help='Allow numbers of mismatch (Default: 3)')

    parser.add_argument('-r', dest='isreversed', action='store_true', default=False, help='crRNA is downstream from PAM (5\'-PAM+gRNA-3\') (Default: False)')

    parser.add_argument('--offind', dest='offind', action='store_true', default=False, help='Do off-target detection using cas-offinder (Default: False)')

    parser.add_argument('--offpath', dest='offpath', type=str, help='cas-offinder program path (with program name)')

    parser.add_argument('--offdir', dest='offdir', type=str, default='', help='Directory only contains the target genome file (Default: same as seqfile)')

    parser.add_argument('--device', dest='device', type=str, choices=['C','G'], default=False, help='Choose CPU or GPU to do cas-offinder (require specific drivers) (Default: auto)')

    parser.add_argument('--threads', dest='threads', type=int, default=1, help='Thread numbers to use, may cause high memory consumption (Default: 1, not use now)')

    parser.add_argument('-o', dest='outfile', type=str, default='results/cas_output.txt', help='Final output file (Default: results/cas_output.txt)')

    return parser


def check_options(parser):

    args = parser.parse_args()

    if not path.exists(path.abspath(path.expanduser(args.seqfile))):

        print("Can not located sequence file, please input full path of sequence file.")

        parser.print_help()

        sys.exit(1)

    if args.offind:

        if not path.exists(path.abspath(path.expanduser(args.offpath))):

            print("Can not located cas-offinder program, please check your path again.")

            parser.print_help()

            sys.exit(1)

    if re.search('[^ACGTRYSWKMBDHVN]', args.pam):

        print("Your PAM is not supported, please change to another PAM sequence.")

        parser.print_help()

        sys.exit(1)

    if not path.exists(path.dirname(path.abspath(path.expanduser(args.outfile)))):

        print('Make directory: %s\n' % path.dirname(path.abspath(path.expanduser(args.outfile))))

        makedirs(path.dirname(path.abspath(path.expanduser(args.outfile))))

    return args


def main():
    
    global work_dir, seqfile, outfile, length, pam, mismatch, isreversed, offdir

    args = check_options(get_options())

    seqfile = path.abspath(path.expanduser(args.seqfile))

    length = args.length

    pam = args.pam

    mismatch = args.mismatch

    isreversed = args.isreversed

    thread = args.threads

    offind = args.offind

    offpath = path.abspath(path.expanduser(args.offpath))

    offdir = args.offdir

    device = args.device

    outfile = args.outfile

    if outfile == sys.stdout:

        work_dir = path.abspath('.')

    else:

        work_dir = path.dirname(path.abspath(path.expanduser(outfile)))

    if offdir == '':

        offdir = path.dirname(seqfile)

    find_putative_seq(seqfile)

    if offind:

        endpoint = cas_offinder(offpath, device)

        if endpoint:

            print('['+datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")+'] ', end=' ')
            print('All finished!\n')


if __name__ == '__main__':

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)


