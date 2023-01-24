import requests

def parse_sequences(file):

    # This is a parsing function for FASTA files

    inFile = open(file, 'r')
    seqList = []
    currentSeq = ''
    for line in inFile:
        if line[0] == ">":
            if currentSeq != '':
                seqList.append(currentSeq)
            currentSeq = ''
        else:
            currentSeq += line.rstrip()
    seqList.append(currentSeq)

    return seqList

def fetch_proteins(URL, Name, num=1):

    # This code is provided by Uniprot URL to fetch proteins --> It gives a FASTA file, which is parsed by the above function

    with requests.get(URL, stream=True) as request:
        request.raise_for_status()
        n = Name + ".fasta"
        with open(n, 'wb') as f:
            for chunk in request.iter_content(chunk_size=2 ** 20):
                f.write(chunk)
    sequencesList = parse_sequences(n)
    if num == 1:
        sequencesList = "".join(sequencesList)
        return sequencesList
    else:
        return sequencesList