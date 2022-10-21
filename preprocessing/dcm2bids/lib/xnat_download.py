#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Created on Fri Aug 27 10:17:31 2021
    Script to download DICOM data from XNAT
    @author: TienTong
"""

import xnat
import argparse
import re

def parseNumList(string):

    m = re.match(r'(\d+)(?:-(\d+))?$', string)
    if not m:
        raise argparse.ArgumentTypeError("'" + string + "' is not a range of number. Expected forms like '0-5' or '2'.")
    start = int(m.group(1))
    if m.group(2): end = int(m.group(2))
    else:          end = start
    return list(range(start, end+1))

def parse_cmdline():

    # Make parser object for cmmand line arguments
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Required arguments
    parser.add_argument("-r", "--range_subject", type=parseNumList,
                        required=True, dest='range_subject',
                        help="specify the subject range to download")
    parser.add_argument("-u", "--username", type=str, required=True,
                        dest='username', help="XNAT username")
    parser.add_argument("-p", "--password", type=str, required=True,
                        dest='password', help="XNAT password")
    parser.add_argument("-d", "--download_dir", type=str, required=True,
                        dest='download_dir', help="Specify download directory")
    return(parser)

def main():

    # Try running with these args, eg: -n 5 -u tongt03
    args = parse_cmdline().parse_args()
    range_subject = args.range_subject
    username = args.username
    password = args.password
    download_dir = args.download_dir

    session  = xnat.connect(server='https://tmii02.mssm.edu/xnat',
                            user=username, password=password, verify=False)
    project  = session.projects["17-1290"]
    subjects = project.subjects

    for i in range_subject:
        subject = subjects[i]
        subject.download_dir(download_dir)

if __name__ == "__main__":
    main()