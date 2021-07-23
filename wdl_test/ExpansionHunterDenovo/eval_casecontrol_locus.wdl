version 1.0

workflow EvalCaseControlLocus {
    input {
        File multisample_profile
        String docker_image
    }
    call RunEHdn {
        input:
            multisample_profile = multisample_profile,
            docker_image = docker_image
    }
    output {
        File data_file = RunEHdn.data_file
    }
}

task RunEHdn {
    input {
        File multisample_profile
        String docker_image
    }
    command {
        ls
    }
    runtime {
        docker: docker_image
    }
    output {
        File data_file = stdout()
    }
}