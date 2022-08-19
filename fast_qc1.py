import os
import gzip, shutil
import math

import numpy as np
import pandas as pd
import matplotlib.patches as patches
import pylab as plt
from Bio import SeqIO





class PreProcessing(object):
    def __init__(self, dir):
        self.dir = dir

    def plot_fastq_qualities(self, filename, ax=None, limit=10000):

        fastq_parser = SeqIO.parse(open(filename, "rt"), "fastq") # to parse fastqfiles
        
        res = [] #empty list to store the quality score of a fastq files
        c = 0
        for record in fastq_parser: # loops through each read
            #print("%s %s" % (record.id, record.seq)) #print each read along with id
            score = record.letter_annotations["phred_quality"] # quality information of 
                                                                #single read in a fastq file
            res.append(score)  #append score of each read into the result
            c += 1
            if c > limit:
                break
        # print("quality score ", res, end="\n")
        df = pd.DataFrame(res) #pandas dataframe is a 2D data array 
                               #like a table with rows & columns..
                               #to load data into DataFrame object(df)
        l = len(df.T) + 1 # length of data in each read
        print("length ", l)

        if ax == None:
            f, ax = plt.subplots(figsize=(12, 5)) # is a function that returns a tuple
            # containing a figure and axes object(s). Thus when using fig, ax = plt.subplots()
            # you unpack this tuple into the variables fig and ax.
        rect = patches.Rectangle((0, 0), l, 20, linewidth=0, facecolor='r', alpha=.4)
                                #create a rectanle patch to a plot
        ax.add_patch(rect) #add the patches to the axes
        rect = patches.Rectangle((0, 20), l, 8, linewidth=0, facecolor='yellow', alpha=.4)
        ax.add_patch(rect) 
        rect = patches.Rectangle((0, 28), l, 12, linewidth=0, facecolor='g', alpha=.4)
        ax.add_patch(rect)
        df.mean().plot(ax=ax, c='blue') #to plot the mean, we have data 
                                        #in the data frame object(df)
        # boxprops = dict(linestyle='-', linewidth=1, color='black')
        # df.plot(kind='box', ax=ax, grid=False, showfliers=False,
        #         color=dict(boxes='yellow', whiskers='black', caps="black"))
        df.boxplot(ax=ax, grid=False, showfliers=False,
                color=dict(boxes='yellow', whiskers='black', caps="black"))


        ax.set_xticks(np.arange(0, l, 5)) #numpy arange() Returns an array with evenly spaced elements as per the interval
        ax.set_xticklabels(np.arange(0, l, 5))

        ax.set_yticks(np.arange(0, 40, 2))
        ax.set_yticklabels(np.arange(0, 40, 2))

        ax.set_xlabel('position')
        ax.set_xlim((0, l))
        ax.set_ylim((0, 40))
        ax.set_title('per base sequence quality')
        ax.set_ylabel('Average Quality Score')
        # ax.set_ylim((0, 200))

        # plt.plot(x_vals, y_vals)

        # plt.ylabel('Average Quality Score')
        # plt.savefig('/Users/neha/Desktop/Data_Files_Python_Analysis/avg_qulaity.png')
        # plt.savefig('/Users/nips/vrisallc/microbe/dummy/zipp/avg_qulity.pdf')
        plt.savefig(filename + ".pdf")
        # plt.show()
        plt.clf()
        return

    def plot_avg_phredscore(self, filename, limit=10000):

        mydict = {}
        c = 0

        for record in SeqIO.parse(open(filename, "rt"), "fastq"):
            quality = record.letter_annotations['phred_quality']

            mydict[c] = quality
            # print('read_' + str(c), mydict[c])
            print('length read_' + str(c), len(mydict[c]))
            while len(mydict[c]) < 251:
                mydict[c].append(0)

            print('length read_' + str(c), len(mydict[c]))
            c += 1
            if c > limit:
                break

        avg_list = []
        sum = 0
        pos_list = []
        # print("dict length ", len(mydict))
        for pos in range(0, 251):
            for j in range(len(mydict)):
                sum += mydict[j][pos]
            # print("sum_", sum)
            avg = sum / len(mydict)
            # print("avg ", avg)
            avg_list.append(avg)
            pos_list.append(pos)
            # avg_list[pos] += avg
            sum = 0

        print("avg list ", avg_list)

        plt.plot(pos_list, avg_list)
        plt.xlabel('Position (bp)')
        plt.ylabel('Average Quality Score')
        # plt.savefig('avg_' + filename + ".pdf")
        plt.savefig('avg1_p1r1.pdf')
        plt.show()
        plt.clf()  #to clear plot variable

    def unzip(self):
        extension = ".gz"
        os.chdir(self.dir)
        for item in os.listdir(self.dir):  # loop through items in dir
            if item.endswith(extension):  # check for ".gz" extension
                gz_name = os.path.abspath(item)  # get full path of files
                # new_file = (os.patchh.basename(gz_name)).rsplit('.', 1)[0]  # get file name for file within
                print('path of input ', gz_name)
                # gz_name = self.files + "/" + item
                new_file = os.path.splitext(item)[0]  #to split the path name into a pair root and ext.
                                                      # 0-gives extension, 1-root path
                # new_file = '/app/outputfile' + '/' + new_file1
                
                # new_file_name = '/Users/nips/vrisallc/python/upload_to_docker/fastqc_vs_code' + '/' + new_file

                with gzip.open(gz_name, "rb") as f_in, open(new_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out) #method to copy the content from source file to destination file

                self.plot_fastq_qualities(new_file, limit=10000)
                self.plot_avg_phredscore(new_file, limit=10000)





#zipped_dir = '/Users/nips/vrisallc/python/fastqc/Input_Files'  # User has to change this with their file location;
# zipped_dir = input("Enter path")
zipped_dir = '/app/files'
pre_process = PreProcessing(zipped_dir)
pre_process.unzip()
