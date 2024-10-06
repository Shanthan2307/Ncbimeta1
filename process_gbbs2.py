import json
import os

import xmltodict
import time
from mdapps import EntrezFetcher





def process_file(input_file_path,output_dir,output_error_dir):
    """Process a single file to extract information."""
    error_list = []
    with open(input_file_path) as json_file:
        json_data = json.load(json_file)

    gb_seqs = json_data.get('GBSet', {}).get('GBSeq', [])
    gb_seqs = gb_seqs if isinstance(gb_seqs, list) else [gb_seqs]
    extracted_infos = []
    num=0
    for gb_seq in gb_seqs:
        
            extracted_infos.append(extract_information(gb_seq,error_list))
            num+=1
    
    # Save extracted data
    base_name = os.path.basename(input_file_path)
    output_file_name = os.path.splitext(base_name)[0] + "_extracted_info.json"
    output_file_path = os.path.join(output_dir, output_file_name)
    if error_list:
        log_errors_to_file(error_list,base_name,output_error_dir)

    with open(output_file_path, 'w') as json_output_file:
        json.dump(extracted_infos, json_output_file, indent=4)

    print(f"Extracted information from {base_name} has been written to {output_file_path}")


def extract_information(gb_seq,error_list):
    """Extract relevant information from the GenBank sequence."""
    first_reference = gb_seq.get("GBSeq_references", {}).get("GBReference", {})
    if isinstance(first_reference, list):
        first_reference = first_reference[0]

    gbseq_xrefs = gb_seq.get('GBSeq_xrefs', {}).get('GBXref', [])
    gbseq_xrefs = gbseq_xrefs if isinstance(gbseq_xrefs, list) else [gbseq_xrefs]

    bs_db_id, bp_db_id = None, None
    for gbxref in gbseq_xrefs:
        if isinstance(gbxref, dict):
            dbname = gbxref.get('GBXref_dbname', '')
            db_id = gbxref.get('GBXref_id', '')
            if dbname == 'BioSample':
                bs_db_id = db_id
            elif dbname == 'BioProject':
                bp_db_id = db_id

    GBSeq_dict = {
        "GBSeq_locus": gb_seq.get("GBSeq_locus", ""),
        "GBSeq_length": gb_seq.get("GBSeq_length", ""),
        "GBSeq_strandedness": gb_seq.get("GBSeq_strandedness", ""),
        "GBSeq_moltype": gb_seq.get("GBSeq_moltype", ""),
        "GBSeq_topology": gb_seq.get("GBSeq_topology", ""),
        "GBSeq_division": gb_seq.get("GBSeq_division", ""),
        "GBSeq_definition": gb_seq.get("GBSeq_definition", ""),
        "GBSeq_source": gb_seq.get("GBSeq_source", ""),
        "GBSeq_organism": gb_seq.get("GBSeq_organism", ""),
        "GBSeq_taxonomy": gb_seq.get("GBSeq_taxonomy", ""),
        "GBReference_title": first_reference.get("GBReference_title", ""),
        "GBReference_journal": first_reference.get("GBReference_journal", ""),
        "GBSeq_comment": gb_seq.get("GBSeq_comment", ""),
        "BioSample_ID": bs_db_id,
        "BioProject_ID": bp_db_id,
        "Gene_classification": 'Hypothetical protein'
    }

    if bs_db_id:
        Biosample_id = EntrezFetcher.search_Biosample(bs_db_id)
        Biosample_details_xml = EntrezFetcher.fetch_Biosample_details(Biosample_id)
        if Biosample_details_xml:
            try:
                data_dict = xmltodict.parse(Biosample_details_xml)
                Biosample_list = extract_additional_info1(data_dict)
                Biosample_dict = {'Biosample' : Biosample_list}
                GBSeq_dict.update(Biosample_dict)
            except Exception as e:
                error_list.append({"biosample_dict": data_dict, "error": str(e)})
                print(f"An error occurred while processing Biosample details: {e}")

    if not GBSeq_dict.get('bioproject_label') and bp_db_id:
        GBSeq_dict['bioproject_label'] = bp_db_id

    return GBSeq_dict

def extract_additional_info1(bio_sample_set):
    """Extract additional information from the BioSample details."""
    # print(bio_sample_set)
    extracted_info = []
    bio_samples = bio_sample_set.get('BioSampleSet', {}).get('BioSample', [])

    # bio_samples = bio_sample_set.get('biosample_dict', {}).get('BioSampleSet', {}).get('BioSample', [])

    # Check if bio_samples is a dictionary; if so, wrap it in a list
    if isinstance(bio_samples, dict):
        bio_samples = [bio_samples]

    for bio_sample in bio_samples:
        info = {}

        # Extract basic fields
        info['taxonomy_name'] = bio_sample.get('Description', {}).get('Organism', {}).get('@taxonomy_name', '')
        info['title'] = bio_sample.get('Description', {}).get('Title', '')
        info['owner_name'] = bio_sample.get('Owner', {}).get('Name', '')
        info['BiosampleId']=bio_sample.get('@accession','')

        # Extract attributes
        attributes = bio_sample.get('Attributes', {}).get('Attribute', [])
        attributes_list = {}
        if isinstance(attributes, dict):
            display_name = attributes.get('@display_name', '')
            text = attributes.get('#text', '')
            if display_name and text:
                attributes_list[display_name] = text
        else:
            for attribute in attributes:
                display_name = attribute.get('@display_name', '')
                text = attribute.get('#text', '')
                if display_name and text:
                    attributes_list[display_name] = text
        info['attributes'] = attributes_list

        # Extract BioProject label
        links = bio_sample.get('Links')
        bioproject_label = ''
        if links:
            links_list = links.get('Link', [])
            if isinstance(links_list, list):
                for link in links_list:
                    if link.get('@target') == 'bioproject':
                        bioproject_label = link.get('@label', '')
                        break
            elif isinstance(links_list, dict) and links_list.get('@target') == 'bioproject':
                bioproject_label = links_list.get('@label', '')
        
        info['bioproject_label'] = bioproject_label
        
        # Add the sample's info to the overall list
        extracted_info.append(info)

    # Output extracted information
    
    return extracted_info

    

def log_errors_to_file(error_inputs, filename,output_error_dir):
    """Log errors to a file for later review."""
    output_file_name = os.path.splitext(filename)[0] + "_error_log.json"
    output_file_path = os.path.join(output_error_dir, output_file_name)

    with open(output_file_path, 'w') as json_output_file:
        json.dump(error_inputs, json_output_file, indent=4)