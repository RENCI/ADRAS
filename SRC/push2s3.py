#!/usr/bin/env python

import os
import logging
from utilities.s3_utilities import utilities as s3_utilities


def main(args):

    if not os.path.exists(args.filename):
        raise IOError(f'Filename {args.filename} does not exist.')

    # the input path (args.path) is the pth for the object in the bucket,
    # not a local path.

    s3_resource = s3_utilities.s3
    logging.debug(s3_resource)
    logging.debug(s3_utilities.config)

    thisBucket=s3_utilities.config['S3_UPLOAD_Main_Bucket']
    thisRegion=s3_utilities.config['region_name']

    logging.info(f'Bucket={thisBucket}')
    logging.info(f'Region={thisRegion}')

    if not s3_utilities.bucket_exists(thisBucket):
        res = s3_utilities.create_bucket(thisRegion, thisBucket)
        logging.info(f'Bucket {thisBucket} created.')
    else:
        logging.info(f'Bucket {thisBucket} already exists.')

    resp = s3_utilities.upload(thisBucket, args.path, args.filename)
    if not resp:
        logging.info(f'Upload to s3://{thisBucket}:/{args.path}/{args.filename} failed.')
    else:
        logging.info(f'Upload to s3://{thisBucket}:/{args.path}/{args.filename} succeeded.')


if __name__ == '__main__':
    from argparse import ArgumentParser
    import sys
    parser = ArgumentParser()
    parser.add_argument('--filename', action='store', dest='filename', type=str, required=True,
                        help='String: Filename to send.')
    parser.add_argument('--path', action='store', dest='path', type=str, required=True,
                        help='String: object path in s3 bucket.')



    args = parser.parse_args()
    sys.exit(main(args))




