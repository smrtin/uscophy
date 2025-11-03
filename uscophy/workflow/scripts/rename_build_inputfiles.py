from pathlib import Path
import sys 
import os 
import shutil


def main(INFILES ,OUT_DIR ):
    counter=0

    try:
        os.mkdir(OUT_DIR)
        print(f"Directory '{OUT_DIR}' created successfully.")
    except FileExistsError:
        print(f"Directory '{OUT_DIR}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{OUT_DIR}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

    for file in INFILES:
        counter+=1 
        filename , extension = os.path.splitext(os.path.basename(file))

        # check if filename matches int - at - int 
        try:
            filename_parts = filename.split("at")
            if len(filename_parts) != 2 :            
                new_filename = str(counter) + 'at100' + extension
                destination = os.path.join(OUT_DIR,new_filename)
                print('{}\t{}'.format(file,destination))
                shutil.copy2(file, destination)
                #print('copy {} to {}'.format(file,destination))
            else:
                
                destination = os.path.join(OUT_DIR,os.path.basename(file))
                shutil.copy2(file, destination)
                #print('copy {} to {}'.format(file,destination))
        except:
            print("An exception occurred") 


if __name__ == "__main__":


    INFILES=list(snakemake.input) 
    OUT_DIR=str(snakemake.output) #output needs to be named in rule!
    main(INFILES,OUT_DIR )
    
    