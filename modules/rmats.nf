process PREPARE_RMATS {
    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    path samples_csv 
    path bams        

    output:
    path "b1.txt", emit: b1
    path "b2.txt", emit: b2
    path "g1_name.txt", emit: g1_name
    path "g2_name.txt", emit: g2_name
    script:
    """
    for bam in *Aligned.sortedByCoord.out.bam; do
        samtools index \$bam
    done
    
    python3 - << 'EOF'
import csv
import glob
import os

groups_dict = {}
with open("${samples_csv}", mode='r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample = row['sample']
        group = row['group']
        if group not in groups_dict:
            groups_dict[group] = []
        groups_dict[group].append(sample)

unique_groups = list(groups_dict.keys())
if len(unique_groups) < 2:
    raise ValueError("The samples.csv file must contain at least two different groups for comparison in rMATS!")

g1_name = unique_groups[0]
g2_name = unique_groups[1]

with open("g1_name.txt", "w") as f_g1:
    f_g1.write(g1_name)
with open("g2_name.txt", "w") as f_g2:
    f_g2.write(g2_name)

print(f"Group 1 (b1): {g1_name} with samples {groups_dict[g1_name]}")
print(f"Group 2 (b2): {g2_name} with samples {groups_dict[g2_name]}")

all_bams = glob.glob("*Aligned.sortedByCoord.out.bam")
bam_map = {}
for b in all_bams:
    for sample in groups_dict[g1_name] + groups_dict[g2_name]:
        if b.startswith(sample):
            bam_map[sample] = os.path.abspath(b)

b1_paths = [bam_map[s] for s in groups_dict[g1_name] if s in bam_map]
b2_paths = [bam_map[s] for s in groups_dict[g2_name] if s in bam_map]

with open("b1.txt", "w") as f1:
    f1.write(",".join(b1_paths))

with open("b2.txt", "w") as f2:
    f2.write(",".join(b2_paths))

print("Generated b1.txt:", ",".join(b1_paths))
print("Generated b2.txt:", ",".join(b2_paths))
EOF
    """
}

process RMATS {
    tag "rMATS_analysis"
    publishDir "${params.outdir}/rmats", mode: 'copy'
    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    path b1_txt
    path b2_txt
    path gtf
    val read_length
    val lib_type

    output:
    path "rmats_output/*", emit: results

    script:
    """
    mkdir -p tmp_rmats rmats_output

    rmats.py \
        --b1 ${b1_txt} \
        --b2 ${b2_txt} \
        --gtf ${gtf} \
        --od rmats_output \
        -t paired \
        --readLength ${read_length} \
        --cstat 0.05 \
        --libType ${lib_type} \
        --tmp tmp_rmats \
        --nthread ${task.cpus}
    """
}

process RMATS2SASHIMIPLOT {
    tag "Sashimi_${event_type}"
    publishDir "${params.outdir}/sashimi_plots/${event_type}", mode: 'copy'
    container 'mordziarz/pipeline_rnaseq:latest' 

    input:
    path b1_txt
    path b2_txt
    path g1_file
    path g2_file
    tuple val(event_type), path(events_txt)

    output:
    path "sashimi_${event_type}_output", emit: plots, optional: true

    script:
    """
    L1=\$(cat ${g1_file} | tr -d '[:space:]')
    L2=\$(cat ${g2_file} | tr -d '[:space:]')
    

    if [ "\$(wc -l < ${events_txt})" -gt 1 ]; then
        
        B1_COUNT=\$(tr ',' '\\n' < ${b1_txt} | grep -v '^\$' | wc -l)
        B2_COUNT=\$(tr ',' '\\n' < ${b2_txt} | grep -v '^\$' | wc -l)
        
        B1_IDX=\$(seq -s, 1 \$B1_COUNT)
        
        B2_START=\$(( B1_COUNT + 1 ))
        B2_END=\$(( B1_COUNT + B2_COUNT ))
        B2_IDX=\$(seq -s, \$B2_START \$B2_END)

        echo "\${L1}: \${B1_IDX}" > groups.txt
        echo "\${L2}: \${B2_IDX}" >> groups.txt

        rmats2sashimiplot \
            --b1 ${b1_txt} \
            --b2 ${b2_txt} \
            -e ${events_txt} \
            --event-type ${event_type} \
            --l1 \${L1} \
            --l2 \${L2} \
            --group-info groups.txt \
            --exon_s 1 \
            --intron_s 5 \
            --font-size 10 \
            --color "#FFA726,#1976D2" \
            -o sashimi_${event_type}_output

    else
        echo "No significant events for ${event_type}, skipping plot generation."
        mkdir -p sashimi_${event_type}_output
    fi
    """
}
