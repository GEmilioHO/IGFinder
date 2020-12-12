
from itertools import islice
import sys
from tqdm import tqdm


# Use the Ensembl REST API: http://rest.ensembl.org/
# We access it from Python with the ensembl_rest library: https://ensemblrest.readthedocs.io/en/latest/
import ensembl_rest

# This is the default client in case the user has not instantiated one
_client = ensembl_rest.EnsemblClient()


# --- HELPER FUNCTIONS

def head(iterable, n=10):
    "Get the first n entries from an iterable."
    iterator = iter(iterable)
    yield from islice(iterator, n)
# ---
        
def chunks_of(iterable, chunk_size=10):
    "Get the entries of the iterable in chunks."
    iterator = iter(iterable)
    while True:
        next_ = list(head(iterator, chunk_size))
        if next_:
            yield next_
        else:
            break
# ---

def divide_region(start, end, max_region_length):
    """ 
    Return a list of subregions of length >= max_region_length.

    The Ensembl API only allows to query regions of 
    length at most 5Mb at a time. 
    """
    subregions = []
    region_upper = start
    while region_upper < end:

        new_region = (region_upper+1, 
                      region_upper + min(max_region_length, end-region_upper))
                      
        subregions.append(new_region)
        
        region_upper += max_region_length

    return subregions
# ---

def unprocessed(arg): return arg
# ---


# --- MAIN LIBRARY FUNCTIONS

def chromosomes_info(species, ensembl_client=_client):
    """
    Fetch the information of the chromosomes of the given species, including it's
    names and lengths.
    """
    # Fetch the genome information using the Ensembl REST API endpoint:
    #    https://rest.ensembl.org/documentation/info/assembly_info
    genome_info = ensembl_client.assembly_info(species=species)
    chrom_data = {chrom['name']: chrom 
                    for chrom in genome_info['top_level_region']
                    if chrom['name'] in genome_info['karyotype']} 

    return chrom_data, genome_info


def overlapping_features(species, region, 
                         feature='gene', 
                         process=unprocessed, 
                         ensembl_client=_client):
    """
    Get the ids of all genes overlapping the given region.

    The process arg receives a function to handle the raw data from the server.
    """
    # Fetch the overlapping features using the Ensembl REST API endpoint:
    #    http://rest.ensembl.org/documentation/info/overlap_region
    raw_data = ensembl_client.overlap_region(species=species, 
                                             region=region, 
                                             params={'feature':feature})
    
    return [process(feature_data) for feature_data in raw_data]
# ---


def genes_in_chrom(species, 
                   chrom, 
                   chrom_length=None, 
                   ensembl_client=_client,
                   max_region_length=5_000_000):
    """
    Fetch info of all the genes in the chromosome.
    """
    if chrom_length is None:
        chrom_length = ensembl_client.assembly_stats(species=species, 
                                                     region=chrom    )['length']
    
    # The Ensembl API only allows to query regions of 
    # length at most 5Mb at a time
    regions_to_query = divide_region(0, chrom_length, max_region_length)
        
    genes = []
    for start, end in tqdm(regions_to_query):
        region = ensembl_rest.region_str(chrom, start, end)

        genes.extend(overlapping_features(species, region, 
                                          process=(lambda gene_info: gene_info['id']),
                                          ensembl_client=ensembl_client))
    
    return genes

# ---


def get_info(feature_ids, 
             process_data=unprocessed, 
             features_per_query=100, 
             ensembl_client=_client):
    """
    From a list of IDs, fetch the available information for each on the server.

    The process_data option receives a function to handle the raw data from 
    the server (default: return as is).
    """
    info = {}
    processed_features = 0
    for chunk in chunks_of(feature_ids, features_per_query):
        # The data must be splitted in chunks to handle the limitations of the 
        # server (max: chunks of 10_000)

        print(f"Progress: {processed_features}/{len(feature_ids)} IDs processed.")

        # Find the data for the chunk of IDs with the help of the Ensembl REST API enpoint:
        #       https://rest.ensembl.org/documentation/info/lookup_post
        new_data = ensembl_client.lookup_post(params={'ids': chunk, 
                                                      'expand': True})

        info.update({gene_id:process_data(gene_data)
            for gene_id,gene_data in new_data.items()
            if gene_data
        })

        processed_features += len(chunk)
        
    return info
# ---
