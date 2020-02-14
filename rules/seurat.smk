rule GO:
    input:
        degFile="results/{contrast}.diffexp.tsv"
    output:
        down="{{contrast}}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt".format(FC = config["FC"],adjp=config["adjp"]),
        up="{{contrast}}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt".format(FC = config["FC"],adjp=config["adjp"])
    params:
        assembly = config["assembly"],
        FC = config["FC"],
        adjp = config["adjp"]
    conda:
        "../envs/runGO.yaml"
    script:
        "../scripts/runGO_singlecell.R"
