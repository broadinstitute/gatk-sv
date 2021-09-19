import os
import numpy as np
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import argparse
# Image helper function
# stack two or more images vertically


def vstack(lst, outname):
    # given a list of image files, stack them vertically then save as
    list_im = lst   # list of image files
    imgs = [Image.open(i) for i in list_im]
    # pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
    min_shape = sorted([(np.sum(i.size), i.size) for i in imgs])[0][1][0]
    imgs_comb = np.vstack([np.asarray(
        i.resize((min_shape, int(i.size[1] / i.size[0] * min_shape)))) for i in imgs])
    # save that beautiful picture
    imgs_comb = Image.fromarray(imgs_comb, "RGB")
    imgs_comb.save(outname)
# combine two images side by side


def hstack(f1, f2, name):
    # given two images, put them side by side, then save to name
    list_im = [f1, f2]
    imgs = [Image.open(i) for i in list_im]
    # pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
    min_shape = sorted([(np.sum(i.size), i.size) for i in imgs])[0][1]
#     print(min_shape)
    imgs_comb = np.hstack([np.asarray(i.resize((min_shape))) for i in imgs])

    # save that beautiful picture
    imgs_comb = Image.fromarray(imgs_comb, "RGB")
    imgs_comb.save(name)


def words(STR1, STR2, outfile, n=100):
    font = ImageFont.truetype("arial.ttf", 70)
    img = Image.new("RGB", (1800, 300), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    draw.text((n, 10), STR1, (0, 0, 0), font=font)
    draw = ImageDraw.Draw(img)
    draw.text((n, 150), STR2, (0, 0, 0), font=font)
    draw = ImageDraw.Draw(img)
    img.save(outfile)

##########
# class Rdplotprefix():
# def __init__(self,variantfile,GetVariantFunc=GetVariants,pedfile,prefixfile,pesrdir,rddir):
# self.variants=GetVariantFunc(inputfile,pedfile,prefixfile).variants


class Variant():
    def __init__(self, chr, start, end, name, type, samples, varname, prefix):
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

    def pesrplotname(self, dir):
        if os.path.isfile(dir + self.varname + ".png"):
            return dir + self.varname + ".png"
        elif os.path.isfile(dir + self.varname + ".left.png") and os.path.isfile(dir + self.varname + ".right.png"):
            hstack(dir + self.varname + ".left.png", dir +
                   self.varname + ".right.png", dir + self.varname + ".png")
            return dir + self.varname + ".png"
        else:
            raise Exception(dir + self.varname + ".png" +
                            " PESR files not found")

    def rdplotname(self, dir, maxcutoff=float("inf")):
        if int(self.end) - int(self.start) > maxcutoff:
            medium = (int(self.end) + int(self.start)) / 2
            newstart = str(round(medium - maxcutoff / 2))
            newend = str(round(medium + maxcutoff / 2))
        else:
            newstart = self.start
            newend = self.end
        if os.path.isfile(dir + self.chr + "_" + newstart + "_" + newend + "_" + self.samples[0] + "_" + self.name + "_" + self.prefix + ".jpg"):
            return dir + self.chr + "_" + newstart + "_" + newend + "_" + self.samples[0] + "_" + self.name + "_" + self.prefix + ".jpg"
        elif os.path.isfile(dir + self.chr + "_" + newstart + "_" + newend + "_" + self.samples[0] + "_" + self.name + "_" + self.prefix + ".jpg"):
            return dir + self.chr + "_" + newstart + "_" + newend + "_" + self.samples[0] + "_" + self.name + "_" + self.prefix + ".jpg"
        else:
            raise Exception(dir + self.chr + "_" + newstart + "_" + newend + "_" +
                            self.samples[0] + "_" + self.name + "_" + self.prefix + ".jpg" + " Rdplot not found")

    def makeplot(self, pedir, rddir, outdir, flank, build="hg38"):
        if self.type != "INS":
            if int(self.end) - int(self.start) < 2000:
                STR2 = self.varname + " " + \
                    str(int(self.end) - int(self.start)) + 'bp'
            else:
                STR2 = self.varname + " " + \
                    str(int((int(self.end) - int(self.start)) / 1000)) + 'kb'
        else:
            STR2 = self.varname
        pesrplot = self.pesrplotname(pedir)
        if self.type == "DUP" or self.type == "DEL":
            rdplot = self.rdplotname(rddir, flank)
            img = Image.open(rdplot)  # rd plot
            # crop out original RD plot annotations
            img2 = img.crop((0, 230, img.size[0], img.size[1]))
            img2.save("croprd.jpg")
            # get new annotation
            STR1 = self.chr + ":" + \
                '{0:,}'.format(int(self.start)) + '-' + \
                '{0:,}'.format(int(self.end)) + " (+" + build + ")"
            outfile = 'info.jpg'
            words(STR1, STR2, outfile, 100)  # new Rd plot
            vstack(['info.jpg', "croprd.jpg", pesrplot], outdir +
                   self.varname + "_denovo.png")  # combine rd pe and sr together
        else:
            STR1 = self.chr + ":" + \
                '{0:,}'.format(int(self.start)) + '-' + \
                '{0:,}'.format(int(self.end)) + " (hg38)"
            outfile = 'info.jpg'
            words(STR1, STR2, outfile, 50)
            vstack(['info.jpg', pesrplot], outdir +
                   self.varname + "_denovo.png")


class VariantInfo():
    def __init__(self, pedfile, prefix):
        self.pedfile = pedfile
        self.prefixdir = {}
        if os.path.isfile(prefix):
            self.prefixfile = prefix
            self.prefix = set([])
            with open(self.prefixfile, "r") as f:
                for line in f:
                    if "#" not in line:
                        prefix, sample = line.rstrip().split()
                        self.prefixdir[sample] = prefix
                        self.prefix.add(prefix)
        else:
            self.prefix = prefix
        famdct = {}
        reversedct = {}
        self.samplelist = []
        with open(pedfile, "r") as f:
            for line in f:
                dat = line.split()
                [fam, sample, father, mother] = dat[0:4]
                if father + "," + mother not in famdct.keys():
                    famdct[father + "," + mother] = [sample]
                else:
                    famdct[father + "," + mother].append(sample)
                reversedct[sample] = father + "," + mother
                self.samplelist.append(sample)
        self.famdct = famdct
        self.reversedct = reversedct
        # QC
        # if self.prefixdir!={}:
        # if set(self.samplelist)!=set(self.prefixdir.keys()):
        # raise Exception("prefix file and ped file has samples mismatch")

    def getprefix(self, sample):
        if self.prefixdir == {}:
            return self.prefix
        else:
            return self.prefixdir[sample]

    def getnuclear(self, sample):
        parents = self.reversedct[sample]
        if parents != "0,0":
            kids = self.famdct[parents].copy()
            kids.remove(sample)
            return sample + ',' + parents
        else:
            return sample


class GetVariants():
    def __init__(self, inputfile, pedfile, prefix):
        self.inputfile = inputfile
        self.variants = []
        self.variantinfo = VariantInfo(pedfile, prefix)
        with open(inputfile, "r") as f:
            for line in f:
                if "#" not in line:
                    dat = line.rstrip().split("\t")
                    [chr, start, end, name, type, samples] = dat[0:6]
                    sample = samples.split(',')[0]
                    varname = samples.split(',')[0] + '_' + name
                    if "," in sample:
                        raise Exception(
                            "should only have 1 sample per variant")
                    prefix = self.variantinfo.getprefix(sample)
                    nuclearfam = self.variantinfo.getnuclear(sample)
                    variant = Variant(chr, start, end, name,
                                      type, nuclearfam, varname, prefix)
                    self.variants.append(variant)

    def GetRdfiles(self):
        with open(self.inputfile + ".igv", "w") as g:
            if self.variantinfo.prefixdir != {}:
                for prefix in self.variantinfo.prefix:
                    open(self.inputfile + '_' + prefix + ".txt", 'w').close()
            else:
                open(self.inputfile + '_' +
                     self.variantinfo.prefix + ".txt", 'w').close()
            for variant in self.variants:
                f = open(self.inputfile + '_' + variant.prefix + ".txt", 'a')
                f.write("\t".join([variant.chr, variant.start, variant.end,
                                   variant.name, variant.type, variant.sample]) + '\n')
                g.write("\t".join([variant.chr, variant.start, variant.end,
                                   variant.name, variant.type, variant.sample, variant.varname]) + '\n')
                f.close()


class GetDenovoPlots():
    def __init__(self, inputfile, pedfile, prefix, pedir, rddir, outdir, flank, build="hg38", GetVariantFunc=GetVariants):
        self.variants = GetVariantFunc(inputfile, pedfile, prefix).variants
        if pedir[-1] == "/":
            self.pedir = pedir
        else:
            self.pedir = pedir + "/"
        if rddir[-1] == "/":
            self.rddir = rddir
        else:
            self.rddir = rddir + "/"
        if outdir[-1] == "/":
            self.outdir = outdir
        else:
            self.outdir = outdir + "/"
        self.build = build
        self.flank = flank

    def getplots(self):
        for variant in self.variants:
            variant.makeplot(self.pedir, self.rddir,
                             self.outdir, self.flank, self.build)

# Main block


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('varfile')
    parser.add_argument('pedfile')
    parser.add_argument('prefix')
    parser.add_argument('flank')
    parser.add_argument('pedir')
    parser.add_argument('rddir')
    parser.add_argument('outdir')
    args = parser.parse_args()
    obj = GetDenovoPlots(args.varfile, args.pedfile, args.prefix, args.pedir,
                         args.rddir, args.outdir, int(args.flank), "hg38", GetVariants)
    obj.getplots()


if __name__ == '__main__':
    main()
