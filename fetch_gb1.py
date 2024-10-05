from Bio import Entrez
import json
import xmltodict
from mpi4py import MPI
import os
from mdapps import EntrezFetcher





def recursive_fetch(query, rank_start, total_articles, chunk_size,mdf,output_error_dir):
    articles_collected = 0
    times = 0
    
    query_dir = query.replace(" ", "_")
    
    while articles_collected < total_articles:
        id_list = EntrezFetcher.search_nucleotide(query, rank_start + articles_collected, chunk_size)
        if len(id_list) == 0:
            print("No articles found.")
            return
        article_details_xml = EntrezFetcher.fetch_nucleotide_details(id_list).decode("utf-8")
        sfile = os.path.join(query_dir, f'art_{rank_start + articles_collected}.xml')
        # Save the XML to a file
        with open(sfile, "w") as xml_file:
            xml_file.write(article_details_xml)
        # Convert XML to JSON
        with open(sfile) as xml_file:
            data_dict = xmltodict.parse(xml_file.read())
            json_data = json.dumps(data_dict)
            # Write the JSON data to a file
            json_filename = os.path.join(query_dir, f"data_gb_{rank_start + articles_collected}_to_{rank_start + articles_collected + len(id_list)}.json")
            with open(json_filename, "w") as json_file:
                json_file.write(json_data)
        times += 1

        # Delete the XML file after conversion
        os.remove(sfile)
        print(f"XML file {sfile} has been deleted.")
        process_file(json_filename,mdf,output_error_dir)


        
        articles_collected += len(id_list)

# Main function to perform the data mining
def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    query = "E.coli and hypothetical protein"
    md = query+' metadata'
    chunk_size = 5
    total_articles = 10 
    rank_start = rank * total_articles
    
    
    
    # Create a directory for the query if it doesn't exist
    query_dir = query.replace(" ", "_")
    mdf = md.replace(" ", "_")
    output_error_dir = 'error_logs'

    if rank == 0:
        if not os.path.exists(query_dir):
            os.makedirs(query_dir)
        if not os.path.exists(mdf):
            os.makedirs(mdf)
        if not os.path.exists(output_error_dir):
            os.makedirs(output_error_dir)
    
    # Synchronize all rank
    # comm.Barrier()
    
    print(f"Rank {rank} searching nucleotide for query: {query}")
    recursive_fetch(query, rank_start, total_articles, chunk_size,mdf,output_error_dir)

if __name__ == "__main__":
    main()

