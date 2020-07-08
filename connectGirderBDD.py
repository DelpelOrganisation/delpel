#!/usr/bin/env python

# import
import girder_client
import click

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-a', '--apikey', default='XXXXXXX', help='apiKey is a "mechanism" to share authentication and rights on folder')
@click.option('-d', '--destinationlocalfolder', default='.', help='Destination Folder')
def syncPatient_click(apikey, destinationlocalfolder):
    """
    \b
    Synchronize delpel patient folder on local system
    apiKey is a "mechanism" to share authentication and rights on folder.
    User should defined through the web interface (inside user information) an apiKey with the corresponding privileges
    """
    syncPatient(apikey, destinationlocalfolder)

def syncPatient(apiKey, destinationLocalFolder):
    # api Rest url of the warehouse
    url = 'http://humanheart-project.creatis.insa-lyon.fr/database/api/v1'

    # Generate the warehouse client
    gc = girder_client.GirderClient(apiUrl=url)

    # Authentication to the warehouse
    gc.authenticate(apiKey=apiKey)

    # download the delpel  collection in local
    gc.downloadResource('5ef3090773e9f0055751d55d', destinationLocalFolder, resourceType='collection', sync=True)

if __name__ == '__main__':
    syncPatient_click()

