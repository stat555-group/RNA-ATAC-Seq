################################################################################
# Download data 
################################################################################

import wget 

def get_encode_ids(file_name, seq_type):
    with open(file_name, 'r') as f:
        lines = f.read().splitlines()
    
    output = []
    for i in range(len(lines)):
        if seq_type in lines[i]:
            output.append(lines[i+1].split(" ")[1])
            output.append(lines[i+2].split(" ")[1])
    return output
    
def download_data(encode_id, output_dir):
    """Download data from encode.
    """
    url = 'https://www.encodeproject.org/files/{}/@@download/{}.tsv'.format(encode_id, encode_id)
    output_file = os.path.join(output_dir, '{}.tsv'.format(encode_id))
    if not os.path.exists(output_file):
        print("Downloading {}".format(url))
        wget.download(url, out=output_file)


if __name__ == '__main__':
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Download data from ENCODE.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file with ENCODE IDs')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output directory')
    parser.add_argument('-t', '--type', type=str, required=True, help='Data type')
    args = parser.parse_args()

    encode_ids = get_encode_ids(args.input, args.type)
    for encode_id in encode_ids:
        download_data(encode_id, args.output)