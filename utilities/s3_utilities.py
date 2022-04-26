#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
#import logging
import pandas as pd

import boto3
from botocore.exceptions import ClientError

class Utilities:

    def __init__(self):
        self.config = self.load_config()
        self.s3 = boto3.resource('s3')

    def load_config(self):
        #config['Access key ID'] = 'XXXXXXXXXXXXXXXXXXXXXXXXX,'
        #config['Secret access key'] = 'XXXXXXXXXXXXXXXXXXXXXXX'
        config = {}
        #f = os.path.join(os.path.expanduser("~"), 'aws_adcirc_credentials.csv')
        #if os.path.exists(f):
        #    df = pd.read_csv(f)
        #    config = df.to_dict()
        #else:
        #   print('Failed to load aws cred file. Terminal.')
        #   sys.exit(1)

        config['S3_UPLOAD_Main_Bucket'] = 'adcirc'
        config['region_name'] = 'us-east-2'
        return config

    def bucket_exists(self, bucket):
        #logging.debug()
        return self.s3.Bucket(bucket) in self.s3.buckets.all()

    def create_bucket(self, region_name, bucket_name):
        bucket_response = self.s3.create_bucket(
                          Bucket=bucket_name,
                          CreateBucketConfiguration={'LocationConstraint': region_name})
        return bucket_response

    def upload(self, thisBucket, path, file):
        try:
        # bucket.upload_file(file, key, ExtraArgs={'ACL':'public-read'})
            resp = self.s3.Object(thisBucket, os.path.join(path, file)).upload_file(Filename=file, ExtraArgs={'ACL':'public-read'})
        except ClientError as e:
            return False
        return True

#    def list_buckets():
#        s3_client = boto3.resource('s3')
#        ans=s3_client.list_buckets()['Buckets']
#        for bucket in ans:
#            print(bucket['Name'])

#############################################################
# instance the utilities class on import so that logging is immediately available.
utilities = Utilities()
