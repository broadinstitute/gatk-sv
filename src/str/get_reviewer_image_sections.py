#!/usr/bin/env python3

"""
Adapted script from "https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/call_non_ref_pathogenic_motifs.py"
"""



import argparse
import ast
import collections
import gzip
import os
import re
#import pandas as pd
from pprint import pformat, pprint
import sys


with open(sys.argv[1], 'rt') as f:
    motif1_svg_tag = f.readline()
    motif1_svg_image_contents = f.read()



def get_reviewer_image_section(motif1_svg_image_contents):
    svg_image_contents_without_defs = motif1_svg_image_contents.replace("</svg>", "").split("</defs>")[-1]
    matches = list(re.finditer("<line[^>]+y1=\"(\d+)\"[^>]+#arrow[^>]+>", svg_image_contents_without_defs, re.DOTALL))
    if len(matches) == 0 or len(matches) >2:
        return "NA", "NA", "NA","NA", "NA", "NA",[],[],"NA","NA","NA","NA","NA","NA"
#For STRs with a songle haplotype:
    if len(matches) == 1:
#T,A,G,C_upper: how many time seqeunce is interrupted starting with nucleotide T,A,G or C.
#This part is not perfectly precise. Depends on how the lines are described in image it will count >T<, >A<,>G< or >C< interruptions.
#The output is not the total number of nucleotides.

        section_upper_contents = svg_image_contents_without_defs
        T_upper=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FCA100\">T<",section_upper_contents, re.DOTALL))
        A_upper=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FF6347\">A<",section_upper_contents, re.DOTALL))
        G_upper=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#2F8734\">G<",section_upper_contents, re.DOTALL))
        C_upper=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#393939\">C<",section_upper_contents, re.DOTALL))
        T_upper_count=len(T_upper)
        A_upper_count=len(A_upper)
        G_upper_count=len(G_upper)
        C_upper_count=len(C_upper)
        Tot_upper_interrrupt=A_upper_count+T_upper_count+G_upper_count+C_upper_count
        resultsT_upper = []
        resultsA_upper = []
        resultsG_upper = []
        resultsC_upper = []
        resultsTot_upper= []
#Below fpart is to extract X coordinates of the orange region (STR region)
# and Y coordinate of the start of the sequence reads (the reference region will be excluded in the bash script)
        upper_match = str(matches[0])
        start_upper_X_orange_match = list(re.finditer("x1=\"(\d+)\"",upper_match, re.DOTALL))
        end_upper_X_orange_match = list(re.finditer("x2=\"(\d+)\"",upper_match, re.DOTALL))
        start_upper_Y_match = list(re.finditer("y1=\"(\d+)\"",upper_match, re.DOTALL))
        start_upper_X_orange=start_upper_X_orange_match[0].group(1)
        end_upper_X_orange=end_upper_X_orange_match[0].group(1)
        start_upper_Y=start_upper_Y_match[0].group(1)
#Below for loop is to extract Y coordinates (lines or reads) of the reads that is found interrupted
        for i in range(int(T_upper_count)):
            ContentT_Yaxis=int(T_upper[i].group(3))
#            if ContentT_Yaxis not in resultsT_upper:
#                 resultsT_upper.append(ContentT_Yaxis)
            if ContentT_Yaxis not in resultsTot_upper:
                 resultsTot_upper.append(ContentT_Yaxis)
        for i in range(int(A_upper_count)):
            ContentA_Yaxis=int(A_upper[i].group(3))
#            if ContentA_Yaxis not in resultsA_upper:
#                 resultsA_upper.append(ContentA_Yaxis)
            if ContentA_Yaxis not in resultsTot_upper:
                 resultsTot_upper.append(ContentA_Yaxis)
        for i in range(int(G_upper_count)):
            ContentG_Yaxis=int(G_upper[i].group(3))
#            if ContentG_Yaxis not in resultsG_upper:
#                 resultsG_upper.append(ContentG_Yaxis)
            if ContentG_Yaxis not in resultsTot_upper:
                 resultsTot_upper.append(ContentG_Yaxis)
        for i in range(int(C_upper_count)):
            ContentC_Yaxis=int(C_upper[i].group(3))
#            if ContentC_Yaxis not in resultsC_upper:
#                 resultsC_upper.append(ContentC_Yaxis)
            if ContentC_Yaxis not in resultsTot_upper:
                 resultsTot_upper.append(ContentC_Yaxis)
        Gap_upper=list(re.finditer("</text><line\sx1=\"(\d+)\".y1=\"(\d+)\"\sx2=\"(\d+)\"\sy2=\"(\d+)\"\sstroke=\"black\"\s/>", section_upper_contents, re.DOTALL))
        Gap_upper_count=len(Gap_upper)
        Gap_bottom_count="NA"
        resultsTot_bottom = []
        Tot_bottom_interrrupt= "NA"
        start_bottom_X_orange = "NA"
        end_bottom_X_orange = "NA"
        start_bottom_Y = "NA"
    if len(matches) == 2:
        section_upper_contents= svg_image_contents_without_defs[matches[0].start():matches[1].start()]
        T_upper=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FCA100\">T<",section_upper_contents, re.DOTALL))
        A_upper=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FF6347\">A<",section_upper_contents, re.DOTALL))
        G_upper=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#2F8734\">G<",section_upper_contents, re.DOTALL))
        C_upper=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#393939\">C<",section_upper_contents, re.DOTALL))
        T_upper_count=len(T_upper)
        A_upper_count=len(A_upper)
        G_upper_count=len(G_upper)
        C_upper_count=len(C_upper)
        Tot_upper_interrrupt=A_upper_count+T_upper_count+G_upper_count+C_upper_count
        resultsT_upper = []
        resultsA_upper = []
        resultsG_upper = []
        resultsC_upper = []
        resultsTot_upper= []
        upper_match = str(matches[0])
        start_upper_X_orange_match = list(re.finditer("x1=\"(\d+)\"",upper_match, re.DOTALL))
        end_upper_X_orange_match = list(re.finditer("x2=\"(\d+)\"",upper_match, re.DOTALL))
        start_upper_Y_match = list(re.finditer("y1=\"(\d+)\"",upper_match, re.DOTALL))
        start_upper_X_orange=start_upper_X_orange_match[0].group(1)
        end_upper_X_orange=end_upper_X_orange_match[0].group(1)
        start_upper_Y=start_upper_Y_match[0].group(1)
        for i in range(int(T_upper_count)):
            ContentT_Yaxis=int(T_upper[i].group(3))
#            if ContentT_Yaxis not in resultsT_upper:
#                 resultsT_upper.append(ContentT_Yaxis)
            if ContentT_Yaxis not in resultsTot_upper:
                 resultsTot_upper.append(ContentT_Yaxis)
        for i in range(int(A_upper_count)):
            ContentA_Yaxis=int(A_upper[i].group(3))
#            if ContentA_Yaxis not in resultsA_upper:
#                 resultsA_upper.append(ContentA_Yaxis)
            if ContentA_Yaxis not in resultsTot_upper:
                 resultsTot_upper.append(ContentA_Yaxis)
        for i in range(int(G_upper_count)):
            ContentG_Yaxis=int(G_upper[i].group(3))
#            if ContentG_Yaxis not in resultsG_upper:
#                 resultsG_upper.append(ContentG_Yaxis)
            if ContentG_Yaxis not in resultsTot_upper:
                 resultsTot_upper.append(ContentG_Yaxis)
        for i in range(int(C_upper_count)):
            ContentC_Yaxis=int(C_upper[i].group(3))
#            if ContentC_Yaxis not in resultsC_upper:
#                 resultsC_upper.append(ContentC_Yaxis)
            if ContentC_Yaxis not in resultsTot_upper:
                 resultsTot_upper.append(ContentC_Yaxis)
        Gap_upper=list(re.finditer("</text><line\sx1=\"(\d+)\".y1=\"(\d+)\"\sx2=\"(\d+)\"\sy2=\"(\d+)\"\sstroke=\"black\"\s/>", section_upper_contents, re.DOTALL))
        Gap_upper_count=len(Gap_upper)
        section_bottom_contents = svg_image_contents_without_defs[matches[1].start():]
        resultsT_bottom = []
        resultsA_bottom = []
        resultsG_bottom = []
        resultsC_bottom = []
        resultsTot_bottom= []
        T_bottom=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FCA100\">T<",section_bottom_contents, re.DOTALL))
        A_bottom=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FF6347\">A<",section_bottom_contents, re.DOTALL))
        G_bottom=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#2F8734\">G<",section_bottom_contents, re.DOTALL))
        C_bottom=list(re.finditer("fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#393939\">C<",section_bottom_contents, re.DOTALL))
        T_bottom_count=len(T_bottom)
        A_bottom_count=len(A_bottom)
        G_bottom_count=len(G_bottom)
        C_bottom_count=len(C_bottom)
        Tot_bottom_interrrupt=A_bottom_count+T_bottom_count+G_bottom_count+C_bottom_count
        bottom_match = str(matches[1])
        start_bottom_X_orange_match = list(re.finditer("x1=\"(\d+)\"",bottom_match, re.DOTALL))
        end_bottom_X_orange_match = list(re.finditer("x2=\"(\d+)\"",bottom_match, re.DOTALL))
        start_bottom_Y_match = list(re.finditer("y1=\"(\d+)\"",bottom_match, re.DOTALL))
        start_bottom_X_orange=start_bottom_X_orange_match[0].group(1)
        end_bottom_X_orange=end_bottom_X_orange_match[0].group(1)
        start_bottom_Y=start_bottom_Y_match[0].group(1)
        for i in range(int(T_bottom_count)):
            ContentT_Yaxis=int(T_bottom[i].group(3))
  #          if ContentT_Yaxis not in resultsT_bottom:
  #               resultsT_bottom.append(ContentT_Yaxis)
            if ContentT_Yaxis not in resultsTot_bottom:
                 resultsTot_bottom.append(ContentT_Yaxis)
        for i in range(int(A_bottom_count)):
            ContentA_Yaxis=int(A_bottom[i].group(3))
  #          if ContentA_Yaxis not in resultsA_bottom:
  #               resultsA_bottom.append(ContentA_Yaxis)
            if ContentA_Yaxis not in resultsTot_bottom:
                 resultsTot_bottom.append(ContentA_Yaxis)
        for i in range(int(G_bottom_count)):
            ContentG_Yaxis=int(G_bottom[i].group(3))
 #           if ContentG_Yaxis not in resultsG_bottom:
 #                resultsG_bottom.append(ContentG_Yaxis)
            if ContentG_Yaxis not in resultsTot_bottom:
                 resultsTot_bottom.append(ContentG_Yaxis)
        for i in range(int(C_bottom_count)):
            ContentC_Yaxis=int(C_bottom[i].group(3))
 #           if ContentC_Yaxis not in resultsC_bottom:
 #                resultsC_bottom.append(ContentC_Yaxis)
            if ContentC_Yaxis not in resultsTot_bottom:
                 resultsTot_bottom.append(ContentC_Yaxis)
        Gap_upper=list(re.finditer("</text><line\sx1=\"(\d+)\".y1=\"(\d+)\"\sx2=\"(\d+)\"\sy2=\"(\d+)\"\sstroke=\"black\"\s/>", section_upper_contents, re.DOTALL))
        Gap_upper_count=len(Gap_upper)
        Gap_bottom=list(re.finditer("</text><line\sx1=\"(\d+)\".y1=\"(\d+)\"\sx2=\"(\d+)\"\sy2=\"(\d+)\"\sstroke=\"black\"\s/>", section_bottom_contents, re.DOTALL))
        Gap_bottom_count=len(Gap_bottom)
    return Tot_upper_interrrupt,len(resultsTot_upper), Gap_upper_count,Tot_bottom_interrrupt,len(resultsTot_bottom),Gap_bottom_count, resultsTot_upper,resultsTot_bottom,start_upper_X_orange,end_upper_X_orange,start_upper_Y,start_bottom_X_orange,end_bottom_X_orange,start_bottom_Y

a,b,c,d,e,f,m,n,x,y,z,k,l,p=get_reviewer_image_section(motif1_svg_image_contents)


file = open(sys.argv[2],"w")
file.write("Times_Upper_Interruptions" + "\t"+ "Count_Upper_reads_with_interruptions" + "\t" +
           "Times_Upper_Gaps" + "\t" + "Times_Bottom_Interruptions" + "\t" 
           "Count_Bottom_reads_with_interruptions" + "\t" + "Times_Bottom_Gaps" +
           "\n" + str(a)+ "\t" + str(b)+ "\t" + str(c)+ "\t" + str(d) +"\t"+
           str(e) +"\t" + str(f) + "\n" )
file.close()
file = open(sys.argv[3],"w")
#file.write("Reads_Upper" + "\t"+ "start_upper_X_orange" + "\t" + "end_upper_X_orange" + "\t" + "start_upper_Y"  + "\t" + " Reads_Bottom" + "\t" + "start_bottom_X_orange" + "\t" + "end_bottom_X_orange" + "\t" + "start_bottom_Y" + "\n" +
#           str(m) +"\t" + str(x) + "\t"  + str(y) +"\t"  +str(z) +"\t" +  str(n) + "\t" + str(k) +"\t" +str(l) +"\t" + str(p) +"\n" )
file.write("start_upper_X_orange" + "\t" + "end_upper_X_orange" + "\t" + "start_upper_Y"  + "\t" + "start_bottom_X_orange" + "\t" + "end_bottom_X_orange" + "\t" + "start_bottom_Y" + "\n" +
           str(x) + "\t"  + str(y) +"\t"  +str(z) +"\t" + str(k) +"\t" +str(l) +"\t" + str(p) +"\n" )
file.close()

#with open(sys.argv[4]) as f1, open(sys.argv[2]) as f2, open(sys.argv[3], "w") as f3:
#    for x, y in zip(f1, f2):
#          f3.write(x.strip() + "\t" + y.strip() + '\n')
