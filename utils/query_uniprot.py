import os
import re
import requests
from requests.adapters import HTTPAdapter, Retry

## Adapted from Uniprot API request generated script

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


# url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=Insulin%20AND%20%28reviewed%3Atrue%29&size=500'
def run_query(outputpath, url):
    progress = 0
    mode = 'w' if not os.path.exists(outputpath) else 'a'
    with open(outputpath, mode) as f:
        for batch, _ in get_batch(url):
            lines = batch.text.splitlines()
            if not progress:
                print(lines[0], file=f)
            for line in lines[1:]:
                print(line, file=f)
            progress += len(lines[1:])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog='Query Uniprot')
    parser.add_argument('outputpath')
    parser.add_argument('url')
    args = parser.parse_args()
    # print(args)
    run_query(args.outputpath,args.url)
