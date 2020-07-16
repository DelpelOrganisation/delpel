# Delpel project

## Database:

Data are available in this [repository](http://humanheart-project.creatis.insa-lyon.fr/database/#collection/5ef3090773e9f0055751d55d)
To synchonize these data with our local folder, use:
````
  python connectGirderBDD.py -a XXXXXXXXXXXXXX -d destinationFolder
  python connectGirderBDD.py -h # For help
````
Replace the XXXXXXXXXXX by your user apiKey (Inside "My Account", you can find API key tab)

## JSON file

You can find the metadata in:
````
  delpel.json
````

## Anonymized patient creation:

You can find the script to generate the anonymized patient creation
````
  python createPatient.py -h
````
You can find the corresponding dictionary used to find the RTStruct names.

To be able to run it, you need https://github.com/OpenSyd/syd_algo
The creation of the research ipp with encryptId algorithm is private
