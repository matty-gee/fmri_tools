#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 20:02:04 2021
@author: TienTong
"""

import argparse
from bids import BIDSLayout
import json
import re

def parse_cmdline():
        # Make parser object
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Required arguments
    parser.add_argument("-s", "--subject", type=str,
                        required=True, dest='subject',
                        help="specify the subject id")
    parser.add_argument("-b", "--bids_path", type=str,
                    required=True, dest='bids_path',
                    help="specify the bids path")
    return(parser)

def modify_json(filename, action="append", new_data=None, delete_fields=None):   
    with open(filename,'r+') as file:
        file_data = json.load(file)
        if action == "append":
            file_data.update(new_data)          
        if action == "delete":
            for key in delete_fields:
                file_data.pop(key, None)
        # sets file's current position at offset.
        file.seek(0)
        # convert back to json.
        json.dump(file_data, file, indent = 4)
        file.truncate()
        
def main():

    args = parse_cmdline().parse_args()
    subject = args.subject 
    bids_path = args.bids_path
    
    layout = BIDSLayout(bids_path)
    
    # add IntendedFor for all fmap
    fmap_json_list = layout.get(subject=subject, extension='json', 
                                datatype='fmap', return_type='filename')
    func_nii_list = layout.get(subject=subject, extension='nii.gz', 
                               datatype='func', return_type='filename', 
                               absolute_paths=False)
    fmap_intend_dict=dict()
    fmap_intend_dict["IntendedFor"] = [('/').join(f.split('/')[1:]) for f in func_nii_list] # dont want the highest directory included
    for fmap_json in fmap_json_list:
        modify_json(new_data=fmap_intend_dict, filename=fmap_json, action="append")
        
    # add TaskName for all func
    for task in layout.get_tasks():
        task_json_list = layout.get(subject=subject, task=task, extension='json', 
                                    datatype='func', return_type='filename')
        taskname_dict=dict()
        taskname_dict["TaskName"] = task
        for task_json in task_json_list:
            modify_json(new_data=taskname_dict, filename=task_json, action="append")
            
    # get list of combined func
    task_combined_json_list=[] 
    ME_json_list = layout.get(subject=subject, task=['rest', 'whyhow'], extension='json', 
                                datatype='func', return_type='filename')
    for task_json in ME_json_list:
        # remove all echo information for the optimally combined func         
        if "echo" not in task_json:
            task_combined_json_list.append(task_json)       
    for task_combined_json in task_combined_json_list:
        # get list of json sidecars that need to be removed      
        with open(task_combined_json,'r+') as file:
            file_data = json.load(file)
            r = re.compile(".*Echo*.")
            exclude_fields = list(filter(r.match, file_data.keys()))
            # remove all Echo information of all combined func 
            modify_json(filename=task_combined_json, action="delete", delete_fields=exclude_fields)
             
if __name__ == "__main__":
    main()