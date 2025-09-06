#!/usr/bin/env python3

import sys
import os
import argparse
import cv2
import numpy as np
from PIL import Image, ImageDraw, ImageFont

def hconcat_resize(img_list, interpolation=cv2.INTER_CUBIC):
    """Horizontally concatenate images with same height"""
    h_min = min(img.shape[0] for img in img_list)
    img_list_resize = [cv2.resize(img, (int(img.shape[1] * h_min / img.shape[0]), h_min), interpolation=interpolation)
                      for img in img_list]
    return cv2.hconcat(img_list_resize)

def vconcat_resize(img_list, interpolation=cv2.INTER_CUBIC):
    """Vertically concatenate images with same width"""
    w_min = min(img.shape[1] for img in img_list)
    img_list_resize = [cv2.resize(img, (w_min, int(img.shape[0] * w_min / img.shape[1])), interpolation=interpolation)
                      for img in img_list]
    return cv2.vconcat(img_list_resize)

def words(STR1, STR2, outfile, n=100):
    """Create annotation image with text"""
    img = Image.new('RGB', (n, n), color=(255, 255, 255))
    draw = ImageDraw.Draw(img)
    # Use default font since arial.ttf may not be available
    try:
        font = ImageFont.load_default()
        draw.text((10, 10), STR1, fill=(0, 0, 0), font=font)
        draw.text((10, 30), STR2, fill=(0, 0, 0), font=font)
    except:
        draw.text((10, 10), STR1, fill=(0, 0, 0))
        draw.text((10, 30), STR2, fill=(0, 0, 0))
    img.save(outfile)

class Variant():
    def __init__(self, chr, start, end, name, type, samples, varname, prefix, family_id):
        self.chr = chr
        self.coord = str(chr) + ":" + str(start) + "-" + str(end)
        self.start = start
        self.end = end
        self.name = name
        self.type = type
        self.prefix = prefix
        self.varname = varname
        self.sample = samples
        self.samples = samples.split(",")
        self.family_id = family_id

    def pesrplotname(self, dir):
        """Find PE/SR plot file"""
        plot_file = os.path.join(dir, f"{self.family_id}_{self.name}.png")
        if os.path.isfile(plot_file):
            return plot_file
        
        # Check for split plots (left/right)
        left_file = os.path.join(dir, f"{self.family_id}_{self.name}.left.png")
        right_file = os.path.join(dir, f"{self.family_id}_{self.name}.right.png")
        
        if os.path.isfile(left_file) and os.path.isfile(right_file):
            left = cv2.imread(left_file)
            right = cv2.imread(right_file)
            horizontal_combined = hconcat_resize([left, right])
            cv2.imwrite(plot_file, horizontal_combined)
            return plot_file
        
        return 'Error'

    def rdplotname(self, dir, maxcutoff=float("inf")):
        """Find RD plot file"""
        if int(self.end) - int(self.start) > maxcutoff:
            medium = (int(self.end) + int(self.start)) / 2
            newstart = str(round(medium - maxcutoff / 2))
            newend = str(round(medium + maxcutoff / 2))
        else:
            newstart = self.start
            newend = self.end
        
        # Try different naming patterns
        patterns = [
            f"{self.chr}_{newstart}_{newend}_{self.samples[0]}_{self.name}_{self.prefix}_{self.samples[0]}.jpg",
            f"{self.chr}_{newstart}_{newend}_{self.samples[0]}_{self.name}_{self.prefix}_{self.family_id}.jpg"
        ]
        
        for pattern in patterns:
            plot_file = os.path.join(dir, pattern)
            if os.path.isfile(plot_file):
                return plot_file
        
        return 'Error'

    def makeplot(self, pedir, rddir, outdir, flank, build="hg38"):
        """Combine IGV and RD plots"""
        if self.type not in ["INS", "snv", "indel", "INS:ME:SVA", "INS:ME:LINE1", "INS:ME:ALU"]:
            if int(self.end) - int(self.start) < 2000:
                STR2 = f"{self.varname} {int(self.end) - int(self.start)}bp"
            else:
                STR2 = f"{self.varname} {int((int(self.end) - int(self.start)) / 1000)}kb"
        else:
            STR2 = self.varname

        # Get plot file names
        pesrplot = self.pesrplotname(pedir)
        rdplot = self.rdplotname(rddir, flank)

        if pesrplot != 'Error' and rdplot != 'Error':
            # Load images
            igv = cv2.imread(pesrplot)
            rd = cv2.imread(rdplot)
            
            if igv is not None and rd is not None:
                # Resize IGV image to remove white space
                y1, x1, h1, w1 = 0, 0, 3000, 800
                if igv.shape[0] > w1 and igv.shape[1] > h1:
                    resized_igv = igv[x1:w1, y1:h1]
                else:
                    resized_igv = igv
                
                # Create annotation
                STR1 = f"{self.chr}:{int(self.start):,}-{int(self.end):,} (+{build})"
                
                # Combine horizontally
                img_h_resize = hconcat_resize([resized_igv, rd])
                
                # Save combined plot
                output_file = os.path.join(outdir, f"{self.varname}.png")
                cv2.imwrite(output_file, img_h_resize)
                
        elif pesrplot != 'Error':
            # Only IGV plot available
            igv = cv2.imread(pesrplot)
            if igv is not None:
                output_file = os.path.join(outdir, f"{self.varname}.png")
                cv2.imwrite(output_file, igv)
                
        elif rdplot != 'Error':
            # Only RD plot available
            rd = cv2.imread(rdplot)
            if rd is not None:
                output_file = os.path.join(outdir, f"{self.varname}.png")
                cv2.imwrite(output_file, rd)

class VariantInfo():
    def __init__(self, pedfile, prefix):
        self.prefix = prefix
        self.reversedctfam = {}
        self.family_samples = {}
        
        if os.path.exists(pedfile):
            with open(pedfile, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= 2:
                        family_id, sample_id = fields[0], fields[1]
                        self.reversedctfam[sample_id] = family_id
                        if family_id not in self.family_samples:
                            self.family_samples[family_id] = []
                        self.family_samples[family_id].append(sample_id)

    def getprefix(self, sample):
        return self.prefix

    def getnuclear(self, sample):
        """Get nuclear family samples for a given sample"""
        family_id = self.reversedctfam.get(sample, sample)
        return ",".join(self.family_samples.get(family_id, [sample]))

class GetVariants():
    def __init__(self, inputfile, pedfile, prefix):
        self.inputfile = inputfile
        self.variants = []
        self.variantinfo = VariantInfo(pedfile, prefix)
        
        with open(inputfile, "r") as f:
            for line in f:
                if "#" not in line and line.strip():
                    dat = line.rstrip().split("\t")
                    if len(dat) >= 6:
                        chr, start, end, name, type, samples = dat[0:6]
                        sample = samples.split(',')[0]
                        varname = f"{name}_{sample}"
                        
                        if "," in sample:
                            raise Exception("should only have 1 sample per variant")
                        
                        prefix = self.variantinfo.getprefix(sample)
                        nuclearfam = self.variantinfo.getnuclear(sample)
                        family_id = self.variantinfo.reversedctfam.get(sample, sample)
                        
                        variant = Variant(chr, start, end, name, type, nuclearfam, varname, prefix, family_id)
                        self.variants.append(variant)

class GetDenovoPlots():
    def __init__(self, inputfile, pedfile, prefix, pedir, rddir, outdir, flank, build="hg38"):
        self.inputfile = inputfile
        self.pedfile = pedfile
        self.prefix = prefix
        self.pedir = pedir
        self.rddir = rddir
        self.outdir = outdir
        self.flank = flank
        self.build = build
        
        os.makedirs(outdir, exist_ok=True)
        self.variants_obj = GetVariants(inputfile, pedfile, prefix)

    def getplots(self):
        """Generate combined plots for all variants"""
        for variant in self.variants_obj.variants:
            variant.makeplot(self.pedir, self.rddir, self.outdir, self.flank, self.build)

def main():
    parser = argparse.ArgumentParser(description="Combine IGV and RD plots")
    parser.add_argument('varfile', help='Variant file')
    parser.add_argument('pedfile', help='Pedigree file')
    parser.add_argument('prefix', help='Prefix for naming')
    parser.add_argument('flank', help='Flank size for large variants')
    parser.add_argument('pedir', help='IGV plots directory')
    parser.add_argument('rddir', help='RD plots directory')
    parser.add_argument('outdir', help='Output directory')
    
    args = parser.parse_args()
    
    obj = GetDenovoPlots(args.varfile, args.pedfile, args.prefix, args.pedir, 
                        args.rddir, args.outdir, int(args.flank), "hg38")
    obj.getplots()

if __name__ == '__main__':
    main() 