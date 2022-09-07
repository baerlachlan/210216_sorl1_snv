rule addUmis:
	input:
		R1 = os.path.join("results", rawData_dir, "fastq", "{SAMPLE}{MERGETAG, .*}" + config["paired_end"]["tags"][0] + config["fastq_ext"]),
		R2 = os.path.join("results", rawData_dir, "fastq", "{SAMPLE}{MERGETAG, .*}" + config["paired_end"]["tags"][1] + config["fastq_ext"]),
		UMI = os.path.join("results", rawData_dir, "fastq", "{SAMPLE}{MERGETAG, .*}" + config["umi"]["add_header"]["tag"] + config["fastq_ext"]),
	output:
		R1 = temp(os.path.join("results", addUmis_dir, "fastq", "{SAMPLE}{MERGETAG, .*}" + config["paired_end"]["tags"][0] + config["fastq_ext"])),
		R2 = temp(os.path.join("results", addUmis_dir, "fastq", "{SAMPLE}{MERGETAG, .*}" + config["paired_end"]["tags"][1] + config["fastq_ext"])),
		UMI1 = temp(os.path.join("results", addUmis_dir, "fastq", "{SAMPLE}{MERGETAG, .*}_I1" + config["fastq_ext"])),
		UMI2 = temp(os.path.join("results", addUmis_dir, "fastq", "{SAMPLE}{MERGETAG, .*}_I2" + config["fastq_ext"])),
		html = os.path.join("results", addUmis_dir, "log", "{SAMPLE}{MERGETAG, .*}.html"),
	conda:
		"../envs/gatk.yml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 4000,
		time = "00-01:00:00",
	shell:
		"""
		fastp \
        	-i {input.R1}  \
        	-I {input.UMI}  \
        	-o {output.R1} \
        	-O {output.UMI1} \
        	--html {output.html} \
			--json /dev/null \
        	--umi --umi_loc=read2 --umi_len=8 \
        	-G -Q -A -L -w 1 -u 100 -n 8 -Y 100

   		fastp \
        	-i {input.R2}  \
        	-I {input.UMI}  \
        	-o {output.R2} \
        	-O {output.UMI2} \
			--html /dev/null \
			--json /dev/null \
        	--umi --umi_loc=read2 --umi_len=8 \
        	-G -Q -A -L -w 1 -u 100 -n 8 -Y 100
		"""