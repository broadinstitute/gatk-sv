version 1.0

import "Structs.wdl"

workflow RDataToTsv {
  input {
    Array[File]  rdata_files            # list of .rData files
    String?      rdata_obj_name         # name of the R object inside each file, e.g. "gene.data.reanno.permu"
    String       r_docker               # docker image with R installed
    RuntimeAttr? runtime_attr_override
  }

  scatter (rdata_file in rdata_files) {
    call ConvertRDataToTsv {
      input:
        rdata_file            = rdata_file,
        rdata_obj_name        = rdata_obj_name,
        r_docker              = r_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  output {
    Array[File] tsv_gz_files = ConvertRDataToTsv.tsv_gz
  }
}


task ConvertRDataToTsv {
  input {
    File         rdata_file
    String?      rdata_obj_name = "gene.data.reanno.permu"
    String       r_docker
    RuntimeAttr? runtime_attr_override
  }

  # Strip .rData suffix and use as output prefix
  String out_prefix = basename(rdata_file, ".rData")
  String out_file   = out_prefix + ".tsv.gz"

  command <<<
    set -euo pipefail

    Rscript - <<'REOF'
    rdata_file     <- "~{rdata_file}"
    rdata_obj_name <- "~{rdata_obj_name}"
    out_file       <- "~{out_file}"

    env <- new.env(parent = emptyenv())
    load(rdata_file, envir = env)

    if (!exists(rdata_obj_name, envir = env)) {
      stop(paste0("Object '", rdata_obj_name, "' not found in ", rdata_file,
                  ". Available objects: ", paste(ls(env), collapse = ", ")))
    }

    dat <- get(rdata_obj_name, envir = env)

    if (!is.data.frame(dat) && !is.matrix(dat)) {
      dat <- as.data.frame(dat)
    }

    con <- gzcon(file(out_file, open = "wb"))
    write.table(dat, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    close(con)

    cat(sprintf("Written %d rows x %d cols -> %s\n", nrow(dat), ncol(dat), out_file))
    REOF
  >>>

  output {
    File tsv_gz = out_file
  }

  RuntimeAttr default_attr = object {
    cpu_cores:         2,
    mem_gb:            16,
    disk_gb:           20 + ceil(size(rdata_file, "GiB") * 6),
    boot_disk_gb:      10,
    preemptible_tries: 3,
    max_retries:       1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu:            select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:         select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks:          "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,     default_attr.boot_disk_gb])
    docker:         r_docker
    preemptible:    select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:     select_first([runtime_attr.max_retries,       default_attr.max_retries])
  }
}
