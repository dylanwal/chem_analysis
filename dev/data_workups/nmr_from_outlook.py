"""

Microsoft Reference for COM objects:
https://learn.microsoft.com/en-us/dotnet/api/microsoft.office.interop.outlook.namespace?view=outlook-pia

https://stackoverflow.com/questions/22813814/clearly-documented-reading-of-emails-functionality-with-python-win32com-outlook
"""

import os
import zipfile
import pathlib

import win32com.client

outlook = win32com.client.Dispatch("Outlook.Application").GetNamespace("MAPI")

# get accounts
folders = outlook.Folders
for i, folder in enumerate(folders):
    print(i, folder.Name)

account = folders[0]  # << change index


# get folder
for i, folder in enumerate(account.Folders):
    print(i, folder.Name)
nmr_data_folder = account.Folders[17]  # << change index


# exit()  # << remove once you set correct index

# get email
emails = [email for email in nmr_data_folder.Items if "DW2-10" in email.Subject]  # << change filter word
# emails = [email for email in nmr_data_folder.Items if "DW2-3" in email.Subject]  # << change filter word


# get attachments
attachments = [email.Attachments[0] for email in emails]  # only takes first attachment in email


# save attachments
path = pathlib.Path(r"C:\Users\nicep\Desktop\DW2_10")  # << change save location

if not os.path.exists(path):
    os.makedirs(path)

print("saving:", len(attachments), "attachments")
for attachment in attachments:
    # save zip folder
    file_location = path / attachment.FileName
    attachment.SaveAsFile(file_location)

    # unzip folder
    with zipfile.ZipFile(file_location, 'r') as zip_ref:
        zip_ref.extractall(path / attachment.FileName.strip(".zip"))

    # delete unzip folder
    os.remove(file_location)

print("done")