#!/usr/bin/env python
#-*- coding: utf-8 -*-
 
"""COnvert clustal output to fasta consensus 

Created at Wed Feb 24 22:40:58 2016 by Kimmo Palin <kpalin@merit.ltdk.helsinki.fi>
"""


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Convert clustal consensus to fasta output")
    
    parser.add_argument("CLUSTAL",
                        help="Input files [default:%(default)s]",
                        nargs="+")
    
    parser.add_argument("-V", "--verbose",default=False,const=True,nargs="?",
                        help="Be more verbose with output [and log to a file] [default:%(default)s]" )

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s:%(funcName)s:%(levelname)s:%(message)s')
        if args.verbose!=True:
            log_file_handler = logging.FileHandler(args.verbose)
            log_file_handler.setFormatter(logging.getLogger().handlers[0].formatter)
            logging.getLogger().addHandler(log_file_handler)

    return args

def clustal2seq(fname):
    from Bio import AlignIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    alignment = AlignIO.read(fname,"clustal")

    for record in alignment:
        if record.id == "CONSENS0":
            break
    
    r = SeqRecord(Seq(str(record.seq).replace("-",""),IUPAC.IUPACAmbiguousDNA),
                  id=fname.replace(".fasta","")+"_CONSENSUS")
    return r
                  
if __name__ == '__main__':
    args=main()

    from Bio import SeqIO

    import sys
    SeqIO.write([clustal2seq(fname)  for fname in args.CLUSTAL],
                    sys.stdout,"fasta")

