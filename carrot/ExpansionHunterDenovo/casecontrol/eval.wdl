version 1.0

workflow EvalCaseControlLocus {
    input {
        File multisample_profile
        File multisample_profile_expected
        String docker_image
    }
    call RunEHdn {
        input:
            multisample_profile = multisample_profile,
            multisample_profile_expected = multisample_profile_expected,
            docker_image = docker_image
    }
    output {
        File data_file = RunEHdn.data_file
    }
}

task RunEHdn {
    input {
        File multisample_profile
        File multisample_profile_expected
        String docker_image
    }
    command <<<
        cmp --silent ~{multisample_profile} ~{multisample_profile_expected} && exit 0 || exit 1
    >>>
    runtime {
        docker: docker_image
    }
    output {
        File data_file = stdout()
    }
}