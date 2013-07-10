import argparse, sys
from files import conversion

def tab_to_fasta(f, seq_col, id_cols):
   for fasta in conversion.tab_to_fasta(f, seq_col, id_cols):
       sys.stdout.write(fasta)
    
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--file', dest='file_to_convert',
                        required=True)
    parser.add_argument('--from', dest='from_fmt', required=True)
    parser.add_argument('--to', dest='to_fmt', required=True)
    parser.add_argument('-sc', '--sequence-column', dest='seq_col',
                        default=None, type=int)
    parser.add_argument('-ic', '--id-cols', dest='id_cols',
                        default=None, nargs='+', type=int)

    args = parser.parse_args()

    decrement = lambda x: x - 1
    
    if(not args.id_cols is None):
        args.id_cols = map(decrement, args.id_cols)

    if(not args.seq_col is None):
        args.seq_col -= 1
        
    print(args.id_cols)
    if(args.from_fmt.startswith('tab')):
        if(args.to_fmt.startswith('fasta')):
            tab_to_fasta(args.file_to_convert, args.seq_col,
                         args.id_cols)


if __name__ == "__main__":
    main()
