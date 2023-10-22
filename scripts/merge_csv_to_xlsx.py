import pandas as pd
import sys
import os

writer = pd.ExcelWriter('default.xlsx') # Arbitrary output name
for csvfilename in sys.argv[1:]:
    df = pd.read_csv(csvfilename)
    df.to_excel(writer,sheet_name=os.path.splitext(csvfilename)[0])
writer.save()