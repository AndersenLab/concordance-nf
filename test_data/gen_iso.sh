# Spike in some diversity in the test set

# ISO 1
zcat BGI1-RET2-AB1-trim-1P.fq.gz > BGI2-RET4-LSJ1-trim-1P--temp.fq
zcat BGI1-RET2-AB1-trim-2P.fq.gz > BGI2-RET4-LSJ1-trim-2P--temp.fq

zcat BGI2-RET5-QX1213-trim-1P.fq.gz | head -n 1000 >> BGI2-RET4-LSJ1-trim-1P--temp.fq
zcat BGI2-RET5-QX1213-trim-2P.fq.gz | head -n 1000 >> BGI2-RET4-LSJ1-trim-2P--temp.fq

gzip BGI2-RET4-LSJ1-trim-1P--temp.fq
gzip BGI2-RET4-LSJ1-trim-2P--temp.fq

# ISO 2
zcat BGI1-RET2-AB1-trim-1P.fq.gz > ECA15-MY2693_S44_L002_1P--temp.fq
zcat BGI1-RET2-AB1-trim-2P.fq.gz > ECA15-MY2693_S44_L002_2P--temp.fq

zcat ECA15-MY2693_S44_L002_1P.fq.gz | head -n 1000 >> ECA15-MY2693_S44_L002_1P--temp.fq
zcat ECA15-MY2693_S44_L002_2P.fq.gz | head -n 1000 >> ECA15-MY2693_S44_L002_1P--temp.fq

gzip ECA15-MY2693_S44_L002_1P--temp.fq
gzip ECA15-MY2693_S44_L002_1P--temp.fq

# ISO 2
zcat NIC1604_RET-S88_S47_L001_1P.fq.gz > JU3167_RET-S66_S60_L001_1P--temp.fq
zcat NIC1604_RET-S88_S47_L001_2P.fq.gz > JU3167_RET-S66_S60_L001_2P--temp.fq

zcat ECA23-JU1511-FCHF7JNBBXX_L5_Sub_CHKPEI85216100518_indexAGGCAGAA_GCATTCG_1P.fq.gz | head -n 1000 >> JU3167_RET-S66_S60_L001_1P--temp.fq
zcat ECA23-JU1511-FCHF7JNBBXX_L5_Sub_CHKPEI85216100518_indexAGGCAGAA_GCATTCG_2P.fq.gz | head -n 1000 >> JU3167_RET-S66_S60_L001_2P--temp.fq

gzip JU3167_RET-S66_S60_L001_1P--temp.fq
gzip JU3167_RET-S66_S60_L001_2P--temp.fq


# Command to rename 
rename -f --subst '--temp' '' *