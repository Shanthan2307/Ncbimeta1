from Bio import Entrez
import json
import xmltodict
import time
import os
from http.client import IncompleteRead
from urllib.error import HTTPError
Entrez.api_key = '4972e450c573d64ce1ed96ab80cbe0e0af09'
Entrez.email = "your_email@example.com"


class EntrezFetcher:
    def __init__(self, email):
        Entrez.email = email
    # Function to search nucleotide for a specific query
    def search_nucleotide(query, retstart, retmax):
        attempt = 0
        max_attempts = 100
        while attempt < max_attempts:
            try:
                handle = Entrez.esearch(db="Nucleotide", term=query, retstart=retstart, retmax=retmax)
                record = Entrez.read(handle)
                handle.close()
                return record["IdList"]
            except HTTPError as e:
                if e.code in [429,400]:
                    print(f"HTTP Error 429: Too Many Requests. Attempt {attempt + 1} failed. Retrying after wait.")
                    attempt += 1
                    time.sleep(60)  # Wait for 60 seconds before retrying
                else:
                    raise e
            except Exception as e:
                print(f"Attempt {attempt + 1} failed: {e}")
                attempt += 1
                time.sleep(5)  # Wait for 5 seconds before retrying
        raise Exception("Failed to search nucleotide after 10 attempts")

    def fetch_nucleotide_details(id_list):
        ids = ",".join(id_list)
        attempt = 0
        max_attempts = 100
        while attempt < max_attempts:
            try:
                handle = Entrez.efetch(db="Nucleotide", id=ids, rettype="gb", retmode="xml")
                records = handle.read()
                handle.close()
                return records
            except IncompleteRead as e:
                print(f"Attempt {attempt + 1} failed: IncompleteRead({e.partial})")
                attempt += 1
                time.sleep(5)  # Wait for 5 seconds before retrying
            except HTTPError as e:
                if e.code in [429,400]:
                    print(f"HTTP Error 429: Too Many Requests. Attempt {attempt + 1} failed. Retrying after wait.")
                    attempt += 1
                    time.sleep(60)  # Wait for 60 seconds before retrying
                else:
                    raise e
            except Exception as e:
                print(f"Attempt {attempt + 1} failed: {e}")
                attempt += 1
                time.sleep(5)  # Wait for 5 seconds before retrying
        raise Exception("Failed to fetch details after 10 attempts")

    def search_Fasta(query):
        """Search for FASTA IDs using the query term."""
        attempt = 0
        max_attempts = 10
        while attempt < max_attempts:
            try:
                handle = Entrez.esearch(db="Nucleotide", term=query, retmax=1)
                record = Entrez.read(handle)
                handle.close()
                return record.get("IdList", [])
            except (IncompleteRead, HTTPError) as e:
                attempt += 1
                print(f"Attempt {attempt} failed: {type(e).__name__} - {e}. Retrying in 5 seconds...")
                time.sleep(5)
            except Exception as e:
                print(f"Attempt {attempt + 1} failed: {e}")
                attempt += 1
                time.sleep(5)
        raise Exception("Failed to search for Biosample after 10 attempts")

    def fetch_Fasta_details(id_list):
        """Fetch details for a list of Fasta IDs."""
        ids = ",".join(id_list)
        attempt = 0
        max_attempts = 10
        while attempt < max_attempts:
            try:
                handle = Entrez.efetch(db="Nucleotide", id=ids, rettype="fasta", retmode="xml")
                records = handle.read()
                handle.close()
                return records
            except (IncompleteRead, HTTPError) as e:
                attempt += 1
                print(f"Attempt {attempt} failed: {type(e).__name__} - {e}. Retrying in 5 seconds...")
                time.sleep(5)
            except Exception as e:
                print(f"Attempt {attempt} failed: {e}")
                attempt += 1
                time.sleep(5)
        raise Exception("Failed to fetch details after 10 attempts")

    def search_Biosample(query):
        """Search for Biosample IDs using the query term."""
        attempt = 0
        max_attempts = 10
        while attempt < max_attempts:
            try:
                handle = Entrez.esearch(db="Biosample", term=query, retmax=100)
                record = Entrez.read(handle)
                handle.close()
                return record.get("IdList", [])
            except (IncompleteRead, HTTPError) as e:
                attempt += 1
                print(f"Attempt {attempt} failed: {type(e).__name__} - {e}. Retrying in 5 seconds...")
                time.sleep(5)
            except Exception as e:
                print(f"Attempt {attempt + 1} failed: {e}")
                attempt += 1
                time.sleep(5)
        raise Exception("Failed to search for Biosample after 10 attempts")

    def fetch_Biosample_details(id_list):
        """Fetch details for a list of Biosample IDs."""
        ids = ",".join(id_list)
        attempt = 0
        max_attempts = 10
        while attempt < max_attempts:
            try:
                handle = Entrez.efetch(db="Biosample", id=ids, rettype="gb", retmode="xml")
                records = handle.read()
                handle.close()
                return records
            except (IncompleteRead, HTTPError) as e:
                attempt += 1
                print(f"Attempt {attempt} failed: {type(e).__name__} - {e}. Retrying in 5 seconds...")
                time.sleep(5)
            except Exception as e:
                print(f"Attempt {attempt} failed: {e}")
                attempt += 1
                time.sleep(5)
        raise Exception("Failed to fetch details after 10 attempts")


# Example usage
if __name__ == "__main__":
    api_key = os.getenv("ENTREZ_API_KEY", "")  # Optionally fetch from environment
    email = os.getenv("ENTREZ_EMAIL", "your_email@example.com")  # Optionally fetch from environment
    fetcher = EntrezFetcher(api_key, email)

    # # Example search for nucleotide
    # try:
    #     nucleotide_ids = fetcher.search_nucleotide("Homo sapiens")
    #     print("Nucleotide IDs:", nucleotide_ids)

    #     # Fetch details for those IDs
    #     if nucleotide_ids:
    #         details = fetcher.fetch_nucleotide_details(nucleotide_ids)
    #         print("Nucleotide Details:", details)
    # except Exception as e:
    #     print("An error occurred:", e)
