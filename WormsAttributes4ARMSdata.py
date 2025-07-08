# This code will read an input file with aphiaIDs (at a minumum) and return 
# a set of attributes for those taxa
# The attributes looked for are 
#"Species importance to society", "IUCN RedList Category","IUCN Criteria","IUCN Year Accessed","HELCOM RedList Category","AMBI ecological group","Environmental position"
# This uses WoRMS REST APIs
# The results are returned in a file with the same number of rows as the input file, with the 
# aphiaID and the attributes listed


import requests
import json
import os
import sys
import csv
import time
import re

# Set defaults
# Note that the input file must have an aphiaID, and you must say which column that is found in 
# All other columns are ignored
# Input filename should be specified here
location = os.path.dirname(os.path.realpath(sys.argv[0]))
infile = location + '\SpeciesList18S.csv'
outfile = location + '\SpeciesListAttributes18S_test.csv'
aphiaIDcol = 1 # column with the aphiaIDs in it

# Functions
# The main API we call
def getAphiaAttributesFromID(aphia_id):
    time.sleep(1.5)
    if aphia_id:
        try:
            species_info = requests.get("https://www.marinespecies.org/rest/AphiaAttributesByAphiaID/"+ aphia_id)
            return species_info
        except Exception as e:
            print(e)
            return None
    else:
        return None    

# For these next two:
# To run the API you need to say which attribute ID you want to look for, and that information 
# is a pain to get, to be honest. I dug it out manually and encoded it here as meastypeids
# The attributes are returned in a formatted JSON file, but we need to dig around in that 
# to find the information we want, which just the attribute value
#
# This function will regex search for "measurementTypeID": {meastypeid} and if it finds it, it will find the 
# closest "measurementValue" and print("measurementValue": {measurementValue})
# Note that attribute 23 can have multiple values in it, but the others only one
meastypeids = [23,1,2,3,131,52,67]
def findallmeasids(json_obj, returnobject=None):
    returnobject = {}
    try:
        for measid in meastypeids:
            dict_return = regexfindjson(json_obj, measid, returnobject)
            #for key value in dict_return: add the key value to the returnobject
            for key, value in dict_return.items():
                returnobject[key] = value
        #go over all keys and values in the returnobject and split the value by ; and then take all the unique elements from that array and merge them back by ; and then add that to the returnobject
        new_returnobject = {}
        for key, value in returnobject.items():
            value_splitted = value.split(";")
            unique_values = set(value_splitted)
            unique_values_joined = ";".join(unique_values)
            new_returnobject[key] = unique_values_joined
        return new_returnobject
    except Exception as e:
        print(e)
        
def regexfindjson(json_obj, measid,returnobject=None):
    if returnobject == None:
        object_to_return = {}
    else:
        object_to_return = returnobject
    try:
        #make regex to find "measurementTypeID": measid,
        regex_to_find = '"measurementTypeID": {measid},'.format(measid=measid)
        #text search the json object for the regex
        does_it_exist = re.search(regex_to_find, json_obj)
        #if it exists, find the closest "measurementValue" and print("measurementValue": {measurementValue})
        if does_it_exist:
            try:
                #get the span of the does_it_exist expression
                spanbegin = does_it_exist.span()[0]
                spanend = does_it_exist.span()[1] + 200
                #get the substring of the json object from spanbegin to spanend
                substring = json_obj[spanbegin:spanend]
                #perform regex search for "measurementValue": "{measurementValue}"" , where {measurementValue} is a string
                regex_measurementValue = '"measurementValue": "(.*?)"'
                measurementvalue = re.search(regex_measurementValue, substring)
                if measurementvalue == None:
                    object_to_return[measid] = "None"             
                if measurementvalue:
                    #get the value of the match by taking the span of the match and getting the substring of the match
                    measurementvalue_match = substring[measurementvalue.span()[0]:measurementvalue.span()[1]]
                    true_meas_value = measurementvalue_match .split(":")[1]
                    #replace the " with nothing
                    true_meas_value = true_meas_value.replace('"', '')
                    #check if there is already a value for this measid in the object_to_return if yes then append the value to the object_to_return
                    if measid in object_to_return:
                        object_to_return[measid] = object_to_return[measid] + ";" + true_meas_value
                    else:
                        object_to_return[measid] = true_meas_value
            except Exception as e:
                print("problem when regex found something")
                pass
            # because it existed and we found a value we can run the function again but exclude the spanend from the json_obj and try and find the next one
            new_json_obj = json_obj[measurementvalue.span()[1]:]
            regexfindjson(new_json_obj, measid=measid,returnobject = object_to_return)
        else:
            if measid not in object_to_return:
                object_to_return[measid] = "None"
        return object_to_return
    except Exception as e:
        print("problem found in regexfindjson")
        print(e)

# Now open the input, run the functions, and write to the output
# 
# 1 Open file
# 2 read first row
# 3 find attribute info (list of actual values or Nones) for the aphiaID 
# 4 add that to the dictionary (if not already present)
# 5 write that row (aphiaID) dictionary to file
# 
# for some bizzare reason, when I do the following - to ask for the length of the list - it 
# completely stops the with on line 137, so I had to do this separately and first 
with open(infile, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    length=len(list(reader))
csvfile.close()
 
taxadict={}
with open(outfile, 'w', newline='') as csvfileout:
    writer = csv.writer(csvfileout, delimiter=',')
    writer.writerow(["AphiaID","Species importance to society", "IUCN RedList Category","IUCN Criteria","IUCN Year Accessed","HELCOM RedList Category","AMBI ecological group","Environmental position"])
    i=0
    with open(infile, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        next(reader)
        for row in reader:
            try:
                aphia_id = row[aphiaIDcol]
                if aphia_id not in taxadict.keys(): 
                    species_info = getAphiaAttributesFromID(aphia_id)
                    if species_info.status_code != 404 or species_info.status_code != 204:
                        raw_json = species_info.json()
                        pretty_json = json.dumps(raw_json, indent=4)
                        regexfound = findallmeasids(pretty_json)
                        taxadict[aphia_id] = [regexfound[meastypeids[0]],regexfound[meastypeids[1]],regexfound[meastypeids[2]],regexfound[meastypeids[3]],regexfound[meastypeids[4]],regexfound[meastypeids[5]],regexfound[meastypeids[6]]]
                    else:
                        taxadict[aphia_id] = ["None","None","None","None","None","None","None"]
            except Exception as e:
                #print("exception:",i,"   ",e)
                taxadict[aphia_id] = ["None","None","None","None","None","None","None"]
            writer.writerow([aphia_id,taxadict[aphia_id][0],taxadict[aphia_id][1],taxadict[aphia_id][2],taxadict[aphia_id][3],taxadict[aphia_id][4],taxadict[aphia_id][5],taxadict[aphia_id][6]])
            #print time left by using (length - i) *1.5 seconds but in minutes, but to avoid loads of prints, only do every 10
            if (i % 50==0): print("time left: " + str((length - i) *1.5/60) + " minutes")
            i+=1
csvfile.close()
csvfileout.close()
 