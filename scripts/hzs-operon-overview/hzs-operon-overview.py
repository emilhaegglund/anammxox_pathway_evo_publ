from Bio import SeqIO
import os
import pandas as pd
import sys


operons = {
    # "MCQ4575098.1": "0 Ca. Bathyanammoxibius GC05_55cm",
    # "MCQ4575097.1": "0 Ca. Bathyanammoxibius GC05_55cm",
    # "MCQ4575096.1": "0 Ca. Bathyanammoxibius GC05_55cm",
    # "MCQ4575095.1": "0 Ca. Bathyanammoxibius GC05_55cm",
    # "MBO1224988.1": "1 Ca. Scalindua sediminis Bin_192",
    # "MBO1224989.1": "1 Ca. Scalindua sediminis Bin_192",
    # "MBO1224990.1": "1 Ca. Scalindua sediminis Bin_192",
    # "MBO1224991.1": "1 Ca. Scalindua sediminis Bin_192",
    # "MBO1224992.1": "1 Ca. Scalindua sediminis Bin_192",
    # "MBO1224993.1": "1 Ca. Scalindua sediminis Bin_192",
    # "MBO1224994.1": "1 Ca. Scalindua sediminis Bin_192",
    "KGMJPFBM_01950": "Ca._Scalindua_erythraensis_AMX13",
    "KGMJPFBM_01951": "Ca._Scalindua_erythraensis_AMX13",
    "KGMJPFBM_01952": "Ca._Scalindua_erythraensis_AMX13",
    "KGMJPFBM_01953": "Ca._Scalindua_erythraensis_AMX13",
    "KGMJPFBM_01954": "Ca._Scalindua_erythraensis_AMX13",
    "KGMJPFBM_01955": "Ca._Scalindua_erythraensis_AMX13",
    "KGMJPFBM_01956": "Ca._Scalindua_erythraensis_AMX13",
    "KGMJPFBM_01957": "Ca._Scalindua_erythraensis_AMX13",
    "KGMJPFBM_01958": "Ca._Scalindua_erythraensis_AMX13",
    "KGMJPFBM_01959": "Ca._Scalindua_erythraensis_AMX13",
    "KGMJPFBM_01960": "Ca._Scalindua_erythraensis_AMX13",
    "GAX62877.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62878.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62879.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62880.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62881.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62882.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62883.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62884.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62885.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62886.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62887.1": "Ca._Scalindua_japonica_husup-a2",
    "GAX62888.1": "Ca._Scalindua_japonica_husup-a2",
    # "GAX62547.1": "2b Ca. Scalindua japonica husup-a2",
    # "GAX62548.1": "2b Ca. Scalindua japonica husup-a2",
    # "GAX62549.1": "2b Ca. Scalindua japonica husup-a2",
    # "GAX62550.1": "2b Ca. Scalindua japonica husup-a2",
    # "GAX62551.1": "2b Ca. Scalindua japonica husup-a2",
    # "GAX62552.1": "2b Ca. Scalindua japonica husup-a2",
    # "GAX62553.1": "2b Ca. Scalindua japonica husup-a2",
    # "QII10642.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_1",
    # "QII10643.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_1",
    # "QII10644.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_1",
    # "QII10645.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_1",
    # "QII10646.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_1",
    # "QII10647.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_1",
    # "QII10648.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_1",
    # "QII10649.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_1",
    # "QII10650.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_1",
    # "QII10651.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_1",
    # "QII12194.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_2",
    # "QII12195.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_2",
    # "QII12196.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_2",
    # "QII12197.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_2",
    # "QII12198.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_2",
    # "QII12199.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_2",
    # "QII12200.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_2",
    # "QII12201.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_2",
    # "QII12202.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_2",
    # "QII12203.1": "Ca._Kuenenia_stuttgartiensis_CSTR1_2",
    # "SOH05194.1": "Ca._Kuenenia_stuttgartiensis_MBR1_1",
    # "SOH05195.1": "Ca._Kuenenia_stuttgartiensis_MBR1_1",
    "SOH05196.1": "Ca._Kuenenia_stuttgartiensis_MBR1_1",
    "SOH05197.1": "Ca._Kuenenia_stuttgartiensis_MBR1_1",
    "SOH05198.1": "Ca._Kuenenia_stuttgartiensis_MBR1_1",
    "SOH05199.1": "Ca._Kuenenia_stuttgartiensis_MBR1_1",
    "SOH05200.1": "Ca._Kuenenia_stuttgartiensis_MBR1_1",
    "SOH05201.1": "Ca._Kuenenia_stuttgartiensis_MBR1_1",
    "SOH05202.1": "Ca._Kuenenia_stuttgartiensis_MBR1_1",
    # "SOH05203.1": "Ca._Kuenenia_stuttgartiensis_MBR1_1",
    "SOH06070.1": "Ca._Kuenenia_stuttgartiensis_MBR1_2",
    "SOH06071.1": "Ca._Kuenenia_stuttgartiensis_MBR1_2",
    "SOH06072.1": "Ca._Kuenenia_stuttgartiensis_MBR1_2",
    "SOH06073.1": "Ca._Kuenenia_stuttgartiensis_MBR1_2",
    "SOH06074.1": "Ca._Kuenenia_stuttgartiensis_MBR1_2",
    "SOH06075.1": "Ca._Kuenenia_stuttgartiensis_MBR1_2",
    "SOH06076.1": "Ca._Kuenenia_stuttgartiensis_MBR1_2",
    "SOH06077.1": "Ca._Kuenenia_stuttgartiensis_MBR1_2",
    "SOH06078.1": "Ca._Kuenenia_stuttgartiensis_MBR1_2",
    # "SOH06079.1": "Ca._Kuenenia_stuttgartiensis_MBR1_2",
    # "GAB63744.1": "Ca._Jettenia_caeni",
    # "GAB63745.1": "Ca._Jettenia_caeni",
    # "GAB63746.1": "Ca._Jettenia_caeni",
    # "GAB63747.1": "Ca._Jettenia_caeni",
    # "GAB63748.1": "Ca._Jettenia_caeni",
    # "GAB63749.1": "Ca._Jettenia_caeni",
    # "GAB63750.1": "Ca._Jettenia_caeni",
    # "GAB63751.1": "Ca._Jettenia_caeni",
    # "GAB63752.1": "Ca._Jettenia_caeni",
    # "GAB63753.1": "Ca._Jettenia_caeni",
    # "GAB63754.1": "Ca._Jettenia_caeni",
    "UJS17933.1": "Ca._Jettenia_AM49_1",
    "UJS17934.1": "Ca._Jettenia_AM49_1",
    "UJS17935.1": "Ca._Jettenia_AM49_1",
    "UJS17936.1": "Ca._Jettenia_AM49_1",
    "UJS17937.1": "Ca._Jettenia_AM49_1",
    "UJS17938.1": "Ca._Jettenia_AM49_1",
    "UJS17939.1": "Ca._Jettenia_AM49_1",
    "UJS17940.1": "Ca._Jettenia_AM49_1",
    "UJS17941.1": "Ca._Jettenia_AM49_1",
    "UJS17942.1": "Ca._Jettenia_AM49_1",
    "UJS18998.1": "Ca._Jettenia_AM49_2",
    "UJS18999.1": "Ca._Jettenia_AM49_2",
    "UJS19000.1": "Ca._Jettenia_AM49_2",
    "UJS19001.1": "Ca._Jettenia_AM49_2",
    "UJS19002.1": "Ca._Jettenia_AM49_2",
    "UJS19003.1": "Ca._Jettenia_AM49_2",
    "UJS19004.1": "Ca._Jettenia_AM49_2",
    "UJS19005.1": "Ca._Jettenia_AM49_2",
    # "UJS19006.1": "Ca._Jettenia_AM49_2",
    # "UJS19007.1": "Ca._Jettenia_AM49_2",
    "GAN34134.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN34135.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN34136.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN34137.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN34138.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN34139.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN34140.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN34141.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN34142.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN34143.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN34144.1": "Ca._Brocadia_sinica_JPN1_1",
    "GAN32116.1": "Ca._Brocadia_sinica_JPN1_2",
    "GAN32117.1": "Ca._Brocadia_sinica_JPN1_2",
    "GAN32118.1": "Ca._Brocadia_sinica_JPN1_2",
    "GAN32119.1": "Ca._Brocadia_sinica_JPN1_2",
    "GAN32120.1": "Ca._Brocadia_sinica_JPN1_2",
    "GAN32121.1": "Ca._Brocadia_sinica_JPN1_2",
    "GAN32122.1": "Ca._Brocadia_sinica_JPN1_2",
    "GAN32123.1": "Ca._Brocadia_sinica_JPN1_2",
    "GAN32124.1": "Ca._Brocadia_sinica_JPN1_2",
    # "GAN32125.1": "Ca._Brocadia_sinica_JPN1_2",
    "BBO18367.1": "Ca._Brocadia_pituitae_317325-1_1",
    "BBO18368.1": "Ca._Brocadia_pituitae_317325-1_1",
    "BBO18369.1": "Ca._Brocadia_pituitae_317325-1_1",
    "BBO18370.1": "Ca._Brocadia_pituitae_317325-1_1",
    "BBO18371.1": "Ca._Brocadia_pituitae_317325-1_1",
    "BBO18372.1": "Ca._Brocadia_pituitae_317325-1_1",
    "BBO18373.1": "Ca._Brocadia_pituitae_317325-1_1",
    "BBO18374.1": "Ca._Brocadia_pituitae_317325-1_1",
    "BBO18375.1": "Ca._Brocadia_pituitae_317325-1_1",
    "BBO18376.1": "Ca._Brocadia_pituitae_317325-1_1",
    "BBO16428.1": "Ca._Brocadia_pituitae_317325-1_2",
    "BBO16429.1": "Ca._Brocadia_pituitae_317325-1_2",
    "BBO16430.1": "Ca._Brocadia_pituitae_317325-1_2",
    "BBO16431.1": "Ca._Brocadia_pituitae_317325-1_2",
    "BBO16432.1": "Ca._Brocadia_pituitae_317325-1_2",
    "BBO16433.1": "Ca._Brocadia_pituitae_317325-1_2",
    "BBO16434.1": "Ca._Brocadia_pituitae_317325-1_2",
    "BBO16435.1": "Ca._Brocadia_pituitae_317325-1_2",
    "BBO16436.1": "Ca._Brocadia_pituitae_317325-1_2",
    "BBO16437.1": "Ca._Brocadia_pituitae_317325-1_2",
    "BBO16438.1": "Ca._Brocadia_pituitae_317325-1_2",
    "UJS19429.1": "Ca._Brocadia_AM9_1",
    "UJS19430.1": "Ca._Brocadia_AM9_1",
    "UJS19431.1": "Ca._Brocadia_AM9_1",
    "UJS19432.1": "Ca._Brocadia_AM9_1",
    "UJS19433.1": "Ca._Brocadia_AM9_1",
    "UJS19434.1": "Ca._Brocadia_AM9_1",
    "UJS19435.1": "Ca._Brocadia_AM9_1",
    "UJS19436.1": "Ca._Brocadia_AM9_1",
    "UJS19437.1": "Ca._Brocadia_AM9_1",
    "UJS19438.1": "Ca._Brocadia_AM9_1",
    "UJS21948.1": "Ca._Brocadia_AM9_2",
    "UJS21949.1": "Ca._Brocadia_AM9_2",
    "UJS21950.1": "Ca._Brocadia_AM9_2",
    "UJS21951.1": "Ca._Brocadia_AM9_2",
    "UJS21952.1": "Ca._Brocadia_AM9_2",
    "UJS21953.1": "Ca._Brocadia_AM9_2",
    "UJS21954.1": "Ca._Brocadia_AM9_2",
    "UJS21955.1": "Ca._Brocadia_AM9_2",
    "UJS21956.1": "Ca._Brocadia_AM9_2",
    "UJS21957.1": "Ca._Brocadia_AM9_2",
    "TWU50255.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU50254.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU50253.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU50252.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU50251.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU50250.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU50249.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU50248.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU50247.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU50246.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU50245.1": "Ca._Brocadia_sapporoensis_B188_1",
    "TWU54191.1": "Ca._Brocadia_sapporoensis_B188_2",
    "TWU54192.1": "Ca._Brocadia_sapporoensis_B188_2",
    "TWU54193.1": "Ca._Brocadia_sapporoensis_B188_2",
    "TWU54194.1": "Ca._Brocadia_sapporoensis_B188_2",
    "TWU54195.1": "Ca._Brocadia_sapporoensis_B188_2",
    "TWU54196.1": "Ca._Brocadia_sapporoensis_B188_2",
    "TWU54197.1": "Ca._Brocadia_sapporoensis_B188_2",
    "TWU54198.1": "Ca._Brocadia_sapporoensis_B188_2",
    "TWU54199.1": "Ca._Brocadia_sapporoensis_B188_2",
    "TWU54200.1": "Ca._Brocadia_sapporoensis_B188_2",
    "TWU54201.1": "Ca._Brocadia_sapporoensis_B188_2",
    # "UJS20456.1": "9c. Ca. Brocadia AM9",
    # "UJS20457.1": "9c. Ca. Brocadia AM9",
    # "UJS20458.1": "9c. Ca. Brocadia AM9",
    # "UJS20459.1": "9c. Ca. Brocadia AM9",
    # "UJS20460.1": "9c. Ca. Brocadia AM9",
    # "QQR65899.1": "A10a. Ca. Brocadia sapporoensis Ega_18-Q3-R5-49",
    # "QQR65900.1": "A10a. Ca. Brocadia sapporoensis Ega_18-Q3-R5-49",
    # "QQR65901.1": "A10a. Ca. Brocadia sapporoensis Ega_18-Q3-R5-49",
    # "QQR65902.1": "A10a. Ca. Brocadia sapporoensis Ega_18-Q3-R5-49",
    # "QQR65903.1": "A10a. Ca. Brocadia sapporoensis Ega_18-Q3-R5-49",
    # "QQR65904.1": "A10a. Ca. Brocadia sapporoensis Ega_18-Q3-R5-49",
    # "QQR65905.1": "A10a. Ca. Brocadia sapporoensis Ega_18-Q3-R5-49",
    # "QQR65906.1": "A10a. Ca. Brocadia sapporoensis Ega_18-Q3-R5-49",
}

orthogroup_path = sys.argv[1]
genbank_dir = sys.argv[2]
orthogroup_out = sys.argv[3]
position_data_out = sys.argv[4]

orthogroup_dict = {}
with open(orthogroup_path) as f:
    for line in f:
        line = line.strip()
        line = line.split()
        for protein_id in line[1:]:
            if protein_id in operons.keys():
                orthogroup_dict[protein_id] = line[0].replace(":", "")


orthogroup_df = pd.DataFrame.from_dict(orthogroup_dict, orient="index")
orthogroup_df["feat_id"] = orthogroup_df.index
orthogroup_df.rename(columns={0: "cluster_id"}, inplace=True)
print(orthogroup_df)
orthogroup_df.to_csv(orthogroup_out, sep="\t", index=False)

operon_data = {
    "seq_id": [],
    "feat_id": [],
    "start": [],
    "end": [],
    "strand": [],
    "product": [],
    "orthogroup": [],
}
for path in os.listdir(genbank_dir):
    gbk_path = os.path.join(genbank_dir, path)
    for record in SeqIO.parse(gbk_path, "genbank"):
        for feat in record.features:
            if feat.type == "CDS":
                if "protein_id" in feat.qualifiers:
                    if feat.qualifiers["protein_id"][0] in operons.keys():
                        feat_id = feat.qualifiers["protein_id"][0]

                        operon_data["seq_id"].append(operons[feat_id])
                        operon_data["feat_id"].append(feat_id)
                        operon_data["start"].append(int(feat.location.start))
                        operon_data["end"].append(int(feat.location.end))
                        if feat.location.strand == 1:
                            operon_data["strand"].append("+")
                        else:
                            operon_data["strand"].append("-")
                        operon_data["product"].append(feat.qualifiers["product"][0])
                        operon_data["orthogroup"].append(orthogroup_dict[feat_id])
                elif "locus_tag" in feat.qualifiers:
                    if feat.qualifiers["locus_tag"][0] in operons.keys():
                        print(feat.qualifiers["locus_tag"][0])
                        feat_id = feat.qualifiers["locus_tag"][0]

                        operon_data["seq_id"].append(operons[feat_id])
                        operon_data["feat_id"].append(feat_id)
                        operon_data["start"].append(int(feat.location.start))
                        operon_data["end"].append(int(feat.location.end))
                        if feat.location.strand == 1:
                            operon_data["strand"].append("+")
                        else:
                            operon_data["strand"].append("-")
                        operon_data["product"].append(feat.qualifiers["product"][0])
                        operon_data["orthogroup"].append(orthogroup_dict[feat_id])

operon_df = pd.DataFrame.from_dict(operon_data)
operon_df.to_csv(position_data_out, sep="\t", index=False)
