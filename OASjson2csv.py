import os, sys, argparse, gzip, json 
from collections import defaultdict
import pandas as pd

text = 'Convert OAS JSON  data to CSV one. For a list of options, type --help or -h.\n\nMetadata: Longitudinal, Chain, Author, Isotype, Age, Size_igblastn, Disease, Link, BSource, BType, Size, Species, Vaccine, Subject.\nBasic data:num_errors, name, seq, v, cdr3, original_name, errors, j, fwh1, fwh2, fwh3, fwh4.\n\nExample usage: python3 OASjson2csv.py ../OAS/json ../OAS/csv -bd original_name seq cdr3 -v'

parser = argparse.ArgumentParser(description=text, formatter_class=argparse.RawTextHelpFormatter)
group = parser.add_mutually_exclusive_group()

parser.add_argument('input', type=str, help='path to input JSON directory')

parser.add_argument('-o', '--output', type=str, help='path to output CSV directory')
group.add_argument('-v', '--verbose', action='store_true', help='display information on used options and data')
group.add_argument('-q', '--quiet', action='store_true', help='hide information about used options and data')
parser.add_argument('-md', '--metadata', nargs='+', type=str, help='select metadata options')
parser.add_argument('-bd', '--basic_data', nargs='+', type=str, help='select basic data options')


args = parser.parse_args()

print(args.metadata)
print(args.basic_data)


meta_v = str()
basic_v = str()


try:
    meta_v = ', '.join(args.metadata)
except TypeError:
    meta_v = args.metadata
try:
    basic_v = ', '.join(args.basic_data)
except TypeError:
    basic_v = args.basic_data
    
if args.quiet:
    print('JSON to CSV converting has started ...')
elif args.verbose:
    print('Metadata: {}\nBasic data: {}'.format(meta_v, basic_v))

    
#Fetch all files in directory and subdirectories.
def get_file_paths(directory):
    for dirpath, _, filenames in os.walk(directory):
        for f in filenames:
            if f.endswith('.json.gz'):
                yield os.path.abspath(os.path.join(dirpath, f))
           
        
# Convert JSON to CSV
def json2csv(meta_d, basic_d, input_file, output_dir=None, verbose=None):    
    meta_keys = []
    basic_keys = []
    
    def prepare_options(options, data):
        if type(options) == list:
            data.extend(options)
        else:
            data.append(options)
        
    prepare_options(meta_d, meta_keys)
    prepare_options(basic_d, basic_keys)
                
    meta_keys = [x for x in meta_keys if x is not None]
    basic_keys = [x for x in basic_keys if x is not None]
    
    to_csv = defaultdict(list)
    to_csv.update((k, []) for k in meta_keys + basic_keys if k is not None)
    
    #The first line are the meta entries.
    meta_line = True
    for line in gzip.open(input_file, 'rb'):
        if meta_line == True:
            metadata = json.loads(line)
            meta_line = False
            for m_key in meta_keys:
                if verbose and meta_keys[0] is not None:
                    print('Metadata: {}'.format(metadata.get(m_key)))
                to_csv[m_key].append(metadata.get(m_key))
            continue
    
        #Parse actual sequence data
        basic_data = json.loads(line)
    
        for b_key in basic_keys:
            if verbose and basic_keys[0] is not None:
                print('Basic data: {}'.format(basic_data.get(b_key)))
            to_csv[b_key].append(basic_data.get(b_key))
    
    data_length = max([len(v) for v in to_csv.values()])
    _ = [to_csv.update({m:to_csv.get(m)*data_length}) for m in meta_keys]
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in to_csv.items()]))
    df.index.name = 'Index'
    
    output_file = input_file.replace('.json.gz', '.csv')
    
    if output_dir:
        output_file = os.path.join(output_dir, output_file.split('/')[-1])
    df.to_csv(output_file)
    

for entry, f in enumerate(list(get_file_paths(args.input))):
    if args.verbose:
        print('Entry #{}'.format(entry))
    json2csv(args.metadata, args.basic_data, f, args.output, verbose=args.verbose)