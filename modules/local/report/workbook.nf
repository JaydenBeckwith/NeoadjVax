// One .xlsx per patient for whoever orders the peptides
// (bin/write_workbook.py).
//
// The TSVs and the LaTeX report are right for a pipeline and wrong for the
// person deciding what to synthesise, who wants to sort, filter and
// comment. This gathers whatever the run produced into one file: the
// ranked candidates per MHC class, the designed long peptides, the HLA-LOH
// calls, the presentation-machinery panel, and a legend that spells out
// every flag.
//
// Optional inputs arrive as [] and are simply absent from the workbook,
// so a DNA-only patient with no LOH file still gets a usable sheet.

process WORKBOOK {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(prioritized), path(long_peptides), path(hla_loh), path(gene_panel)

    output:
    tuple val(meta), path("${meta.id}.xlsx"), emit: xlsx

    script:
    def max_rows = params.workbook_max_rows ? "--max-rows ${params.workbook_max_rows}" : ''
    """
    # Sheets are built in the order the wet lab reads them: what to make,
    # then the peptides themselves, then the context behind the calls.
    args=()

    for report in ${prioritized}; do
        case "\$report" in
            *MHC_Class_I.*|*MHC_I.*) name="Class I candidates" ;;
            *MHC_Class_II.*|*MHC_II.*) name="Class II candidates" ;;
            *Combined*) name="All candidates" ;;
            *) name="Candidates" ;;
        esac
        args+=(--sheet "\$name=\$report")
    done

    for peptides in ${long_peptides}; do
        [ -e "\$peptides" ] || continue
        case "\$peptides" in
            *MHC_Class_I.*|*MHC_I.*) name="Long peptides class I" ;;
            *MHC_Class_II.*|*MHC_II.*) name="Long peptides class II" ;;
            *) name="Long peptides" ;;
        esac
        args+=(--sheet "\$name=\$peptides")
    done

    for loh in ${hla_loh}; do
        [ -e "\$loh" ] || continue
        args+=(--sheet "HLA LOH=\$loh")
    done

    for panel in ${gene_panel}; do
        [ -e "\$panel" ] || continue
        args+=(--sheet "Gene panel=\$panel")
    done

    write_workbook.py "\${args[@]}" ${max_rows} -o ${meta.id}.xlsx
    """

    stub:
    """
    # A real .xlsx is a zip, so the stub writes an empty one rather than
    # a text file pretending to be a workbook.
    python3 -c "import zipfile; zipfile.ZipFile('${meta.id}.xlsx', 'w').close()"
    """
}
