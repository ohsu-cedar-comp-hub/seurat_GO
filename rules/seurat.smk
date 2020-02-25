rule GO:
    input:
        degFile="input/{contrast}.diffexp.tsv"
    output:
        down="results/{contrast}/{contrast}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt",
        up="results/{contrast}/{contrast}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt"
    params:
        assembly = config["assembly"],
        FC = config["FC"],
        adjp = config["adjp"],
        out_dir = lambda w: "results/{}".format(w.contrast),
        up_barplot_out = lambda w: "results/{}/{}.upFC.{}.adjp.{}.BP_GO_barplot.pdf".format(w.contrast, w.contrast, w.FC, w.adjp),
        up_dag_out = lambda w: "results/{}/{}.upFC.{}.adjp.{}.BP_GO_dag".format(w.contrast, w.contrast, w.FC, w.adjp),
        down_barplot_out = lambda w: "results/{}/{}.downFC.{}.adjp.{}.BP_GO_barplot.pdf".format(w.contrast, w.contrast, w.FC, w.adjp),
        down_dag_out = lambda w: "results/{}/{}.downFC.{}.adjp.{}.BP_GO_barplot.pdf".format(w.contrast, w.contrast, w.FC, w.adjp),
        up_consolidated_out = lambda w: "results/{}/{}.upFC.{}.adjp.{}.BP_GO_consolidated.tsv".format(w.contrast, w.contrast, w.FC, w.adjp),
        down_consolidated_out = lambda w: "results/{}/{}.downFC.{}.adjp.{}.BP_GO_consolidated.tsv".format(w.contrast, w.contrast, w.FC, w.adjp),
    conda:
        "../envs/runGO.yaml"
    script:
        "../scripts/runGO_singlecell.R"
