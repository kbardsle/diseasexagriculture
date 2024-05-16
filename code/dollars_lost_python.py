#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:38:35 2024

@author: emmakeppler
"""

# Python  script for processing crop total production data from NASS quickstats
# Download NASS crop total data prior to running script

import pandas as pd

# open csv in pandas
crop_dollars = pd.read_csv('crop_production_dollars.csv')

# only keep state, value, category columns
# only take $ sales category
crop_dollars = crop_dollars.loc[crop_dollars['Data Item']=='CROP TOTALS - SALES, MEASURED IN $',
            ['State', 'Domain Category','Value']]

# NAICS classification starting with 111 indicates crops not other agricultural products
crop_dollars = crop_dollars[crop_dollars['Domain Category'].str.contains('NAICS CLASSIFICATION: (111', regex=False)] 

# remove non-number values (missing data), strip commas, and convert number strings to float
crop_dollars = crop_dollars[~crop_dollars['Value'].str.contains('(', regex=False)]
crop_dollars['Value']= crop_dollars['Value'].str.replace(',','').astype('float')

# sum all row values for sales by state into new dataframe
dollar_list = []
for s in crop_dollars['State'].unique():
    dollar = crop_dollars.loc[crop_dollars['State']==s,'Value'].sum()
    dollar_list.append([s,dollar])
state_production = pd.DataFrame(data=dollar_list,  columns = ['state','dollars'])


# load hours lost csv (run hours lost script)
hours = pd.read_csv('state_hours_lost.csv')

# match keys for merge
hours['State'] = hours['State'].str.lower()
state_production['state'] = state_production['state'].str.lower()
hours = hours.rename({'State':'state'}, axis=1)

# merge dataframes on state
pop_dollars = pd.merge(hours, state_production, on='state')

# multiply proportion hours lost by dollars to get lost dollars
pop_dollars.insert(len(pop_dollars.columns), column = 'lost_dollars', 
                   value = pop_dollars['prop_yearly_hrs_lost']*pop_dollars['dollars'])

# save merged dataframe with lost dollars
pop_dollars.to_csv('state_lost_dollars.csv', index=False)
