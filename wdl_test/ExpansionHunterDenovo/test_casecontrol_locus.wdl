version 1.0

workflow TestCaseControlLocus {
    input {
        Int test
        String docker_image
    }
    call RunEHdn {
        input:
            test = test,
            docker_image = docker_image
    }
    output {
        File data_file = RunEHdn.data_file
    }
}

task RunEHdn {
    input {
        Int test
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