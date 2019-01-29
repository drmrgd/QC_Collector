#!/bin/bash
# Super quickie script to uncompress and unarchive IR data when downloaded directly from the 
# IR REST API.  Took way longer to finally get around to writing this than it should!
##################################################################################################
VERSION='1.5.0_062317'
cwd=$(pwd)

required_progs=('parallel' 'rename')
for p in ${required_progs[@]}; do
    command -v $p > /dev/null 2>&1 || {
        echo "ERROR: $p is not found on this system but is required. Please install '$p' before proceeding."
        exit 1
    }
done

# Check for download_zips dir and make one if need be
echo -n "Checking for a 'download_zips' directory..."
if [[ ! -e $cwd/download_zips ]]; then
    echo -e "\n\tNo 'download_zips' directory found. Creating new."
    mkdir "$cwd/download_zips"
else
    echo "Done!"
fi

# Use GNU Parallel to open the packages
echo -n "Unzipping IR archive data..."
if [[ ! $(find . -maxdepth 1 -iname "*download.zip") ]]; then 
    echo -e "\nERROR: No IR API *download.zip file(s) can be found in this directory!"
    exit 1
else
    parallel unzip -q {} > /dev/null 2>&1 ::: "$cwd/*zip"
    if [[ $? -ne 0 ]]; then
        echo -e "\nERROR: There was a problem unzipping the IR archive files!"
        exit 1
    else
        echo "Done!"
    fi
fi

# Move the download.zip files to the 'download_zips' dir
echo -n "Moving download.zip files to archive directory, and removing log files..."
mv *download.zip $cwd/download_zips/
rm *log
echo "Done!"

# Get the sample name and unzip each into its own new directory
echo "Unzipping IR data..."
for zip in *zip; do
    name=''
    if [[ $(echo $zip | egrep '.*_c[0-9]{3,}_.*') ]]; then
        name=$(echo $zip | perl -pe 's/_c[0-9]{3,}_.*//')
    elif [[ $(echo $zip | egrep '_[0-9a-f]{8}-[0-9a-f]{4}.*') ]]; then
        name=$(echo $zip | perl -pe 's/^(.*?(?:_v[0-9]+)?)_[0-9a-f]{8}-.*/$1/')
    else
        #Need a fallback to handle all of these odd cases
        name=$(echo $zip | sed 's/\.zip//')
    fi
    echo -e "\t Processing $name..."
    unzip -q -d $name $zip && rm $zip
done
echo "Done!"

# Create a VCFs directory
echo -n "Checking for a 'vcfs' directory to collect the VCF files..."
if [[ ! -e $cwd/vcfs ]]; then
    echo -e "\n\tNo 'vcfs' directory found. Creating a new one."
    mkdir "$cwd/vcfs"
elif [[ $(find "$cwd/vcfs" -maxdepth 1 -name "*vcf") ]]; then
    echo -e "\n\tWARN: There are VCF files already in this directory which may be overwritten by newer version!"
    echo "Done!"
fi

# Copy the VCF files into the new directory
echo -n "Copying VCF files into 'vcfs' directory and trimming name..."
find . -iname "*vcf" -not -name 'SmallVariants*vcf' -not -name '*_Filtered_*' -exec cp {} "$cwd/vcfs/" \; > /dev/null 2>&1

# Fix stupid name from IR5.2. If there is VCF here from an older analyis already when we do this, overwrite it. We always
# want the latest file anyway.  
cd vcfs
for vcf in *vcf; do
    if [[ $(echo $vcf | egrep '_Non-Filtered_201[78]-[0-9]{2}.*vcf') ]]; then
        new_name=${vcf/_Non-Filtered*/.vcf}
        mv $vcf $new_name
    fi
done
echo "Done!"
