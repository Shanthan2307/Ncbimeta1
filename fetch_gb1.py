import multiprocessing as mp
from Bio import Entrez
import json
import xmltodict
import os
import time  # Import the time module to measure execution time
from mdapps import EntrezFetcher
from process_gbbs2 import process_file

# Define the recursive fetch function
def recursive_fetch(query, rank_start, total_articles, chunk_size, mdf, output_error_dir):
    articles_collected = 0
    query_dir = query.replace(" ", "_")
    
    while articles_collected < total_articles:
        id_list = EntrezFetcher.search_nucleotide(query, rank_start + articles_collected, chunk_size)
        if len(id_list) == 0:
            print(f"No articles found for rank {rank_start}.")
            return
        article_details_xml = EntrezFetcher.fetch_nucleotide_details(id_list).decode("utf-8")
        
        # Save the XML to a file
        sfile = os.path.join(query_dir, f'art_{rank_start + articles_collected}.xml')
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
        
        # Delete the XML file after conversion
        os.remove(sfile)
        print(f"XML file {sfile} has been deleted.")
        
        # Process the JSON file
        process_file(json_filename, mdf, output_error_dir)
        
        articles_collected += len(id_list)

# Main function to manage multiprocessing
def main():
    start_time = time.time()  # Start timer for total execution time
    
    query = "E.coli and hypothetical protein"
    md = query + ' metadata'
    chunk_size = 50
    total_articles = 200
    num_processes = 10  # Define the number of processes you want to run
    query_dir = query.replace(" ", "_")
    mdf = md.replace(" ", "_")
    output_error_dir = 'error_logs'

    # Create directories if they don't exist
    if not os.path.exists(query_dir):
        os.makedirs(query_dir)
    if not os.path.exists(mdf):
        os.makedirs(mdf)
    if not os.path.exists(output_error_dir):
        os.makedirs(output_error_dir)
    
    # Define a pool of workers for multiprocessing
    with mp.Pool(processes=num_processes) as pool:
        tasks = []
        
        # Distribute the work across processes
        for rank_start in range(0, total_articles * num_processes, total_articles):
            tasks.append(pool.apply_async(recursive_fetch, (query, rank_start, total_articles, chunk_size, mdf, output_error_dir)))
        
        # Collect results from the workers
        for task in tasks:
            task.get()  # Ensure each task completes
    
    end_time = time.time()  # End timer for total execution time
    total_time_seconds = end_time - start_time
    total_time_minutes = total_time_seconds / 60  # Convert time to minutes
    print(f"Total time taken to complete the script: {total_time_minutes:.2f} minutes")

if __name__ == "__main__":
    main()
