import os
import pytest

from sv_utils import fix_vcf, genomics_io


class Default:
    resources_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
    small_needs_fix_vcf = os.path.join(resources_dir, "small_needs_fix.vcf.gz")


def gatk_installed() -> bool:
    return os.system("gatk") == 0


def vcf_is_valid(vcf: str) -> bool:
    return os.system(f"gatk ValidateVariants -V {vcf} > /dev/null 2>&1") == 0


@pytest.mark.integration_test
def test_fix_vcf(
        tmpdir,
        vcf: str = Default.small_needs_fix_vcf
):
    # The whole point of fix_vcf is to make the VCF readable by GATK. However, I don't want the test to mean that GATK
    # must be installed on any system to pass tests. So if GATK is not installed, just run through the code to be sure
    # that it
    if gatk_installed():
        # original VCF should not be valid
        assert not vcf_is_valid(vcf)
    # fix the vcf
    temp_out_dir = tmpdir.mkdir("test_fix_vcf")
    fixed_vcf = os.path.join(temp_out_dir, os.path.basename(vcf))
    fix_vcf.fix_vcf(input_vcf=vcf, output_vcf=fixed_vcf)
    if gatk_installed():
        # now it should be valid
        assert vcf_is_valid(fixed_vcf)
    else:
        # make sure that the file exists at that path, and is readable as a VCF
        assert os.path.isfile(fixed_vcf)
        genomics_io.vcf_to_pandas(fixed_vcf)
