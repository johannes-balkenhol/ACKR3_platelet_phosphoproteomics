###################################################################################################
###################################################################################################
########### import and process top.table.pka.mita with specified phospho site for flagging early and late timepoints
use lipidomics;
### use import wizared in mysql workbench for tables 
### top.all.csv, top.pka.mita.csv
### specified_phosphosite.csv, norm_abundance.csv, raw_abundance.csv

####################################################################################################
####################################################################################################
### annotation of data and filtering


### annotate specified phosphosites
alter table `lipidomics`.`top.pka.mita` add column specified_phosphosite int(2);
Update `lipidomics`.`top.pka.mita` set specified_phosphosite = 0;

update `lipidomics`.`top.pka.mita` a
inner join `lipidomics`.`specified_phosphosite` b
on a.ID = b.ID
set a.specified_phosphosite = 1;



### annotate description
alter table `lipidomics`.`top.pka.mita` add column description mediumtext;

update `lipidomics`.`top.pka.mita` a
inner join phosphop b
on a.uniprot_id_10 = b.Uniprot_id
set a.description = b.Description;
### not getting all descriptions


### annote the timepoint regulation
alter table `lipidomics`.`top.pka.mita` add column early_response int(2);
alter table `lipidomics`.`top.pka.mita` add column late_response int(2);

update `lipidomics`.`top.pka.mita` set early_response = 0;
update `lipidomics`.`top.pka.mita` set late_response = 0;

update `lipidomics`.`top.pka.mita` set early_response = 1 
where `adj.p.value_0010` <= 0.05
or `adj.p.value_0030` <= 0.05
or `adj.p.value_0060` <= 0.05
or `adj.p.value_0300` <= 0.05;

update `lipidomics`.`top.pka.mita` set late_response = 1 
where `adj.p.value_0600` <= 0.05
or `adj.p.value_0900` <= 0.05
or `adj.p.value_1800` <= 0.05;



### check or different hypotheis
## different response times
## all with residue
SELECT distinct ID,uniprot_id_10 as uniprot_id,description,symbol_0010 as symbol, Residue,Site, Sequence, specified_phosphosite,
logfc_0010,logfc_0030,logfc_0060,logfc_0300,logfc_0600,logfc_0900,logfc_1800, early_response,late_response,
`adj.p.value_0010`,`adj.p.value_0030`,`adj.p.value_0060`,`adj.p.value_0300`,`adj.p.value_0600`,`adj.p.value_0900`,`adj.p.value_1800`
FROM lipidomics.`top.pka.mita`
where Residue != ""
order by symbol asc;

### group by the identified sites
SELECT distinct ID,uniprot_id_10 as uniprot_id,description,symbol_0010 as symbol, Residue,Site, Sequence, specified_phosphosite,
logfc_0010,logfc_0030,logfc_0060,logfc_0300,logfc_0600,logfc_0900,logfc_1800, early_response,middle_response,late_response,
`adj.p.value_0010`,`adj.p.value_0030`,`adj.p.value_0060`,`adj.p.value_0300`,`adj.p.value_0600`,`adj.p.value_0900`,`adj.p.value_1800`
FROM lipidomics.`top.pka.mita`
where  late_response = 1
#and specified_phosphosite = 1
and Residue != ""
group by uniprot_id, symbol, residue, site
order by symbol asc;


####################################################################################################
####################################################################################################
######## Resultstble Table 1 
######## signficant Diff. regulated specified phosphosites that are PKA targets
######## uiprot_id, symbol, residue, site, sequence, logfc, time
drop table table1; 
create table table1
(SELECT distinct uniprot_id_10 as uniprot_id,symbol_0010 as symbol, Residue,Site, Sequence, specified_phosphosite, avg(logfc_0010), "10s" as time
FROM lipidomics.`top.pka.mita`
where  `adj.p.value_0010`< 0.05 and Residue != ""
group by uniprot_id, symbol, residue, site order by avg(logfc_0010) desc)
UNION
(SELECT distinct uniprot_id_30 as uniprot_id,symbol_0030 as symbol, Residue,Site, Sequence, specified_phosphosite, avg(logfc_0030), "30s" 
FROM lipidomics.`top.pka.mita`
where  `adj.p.value_0030`< 0.05 and Residue != ""
group by uniprot_id, symbol, residue, site order by avg(logfc_0030) desc)
UNION
(SELECT distinct uniprot_id_30 as uniprot_id,symbol_0060 as symbol, Residue,Site, Sequence, specified_phosphosite, avg(logfc_0060), "60s" 
FROM lipidomics.`top.pka.mita`
where  `adj.p.value_0060`< 0.05 and Residue != ""
group by uniprot_id, symbol, residue, site order by avg(logfc_0060) desc)
UNION
(SELECT distinct uniprot_id_300 as uniprot_id,symbol_0300 as symbol, Residue,Site, Sequence, specified_phosphosite, avg(logfc_0300), "300s" 
FROM lipidomics.`top.pka.mita`
where  `adj.p.value_0300`< 0.05 and Residue != ""
group by uniprot_id, symbol, residue, site order by avg(logfc_0300) desc)
UNION
(SELECT distinct uniprot_id_600 as uniprot_id,symbol_0600 as symbol, Residue,Site, Sequence, specified_phosphosite, avg(logfc_0600), "600s" 
FROM lipidomics.`top.pka.mita`
where  `adj.p.value_0600`< 0.05 and Residue != ""
group by uniprot_id, symbol, residue, site order by avg(logfc_0600) desc)
UNION
(SELECT distinct uniprot_id_900 as uniprot_id,symbol_0900 as symbol, Residue,Site, Sequence, specified_phosphosite, avg(logfc_0900), "900s" 
FROM lipidomics.`top.pka.mita`
where  `adj.p.value_0900`< 0.05 and Residue != ""
group by uniprot_id, symbol, residue, site order by avg(logfc_0900) desc)
UNION ALL
(SELECT distinct uniprot_id_1800 as uniprot_id,symbol_1800 as symbol, Residue,Site, Sequence, specified_phosphosite, avg(logfc_1800), "1800s" 
FROM lipidomics.`top.pka.mita`
where  `adj.p.value_1800`< 0.05 and Residue != ""
group by uniprot_id, symbol, residue, site order by avg(logfc_1800) desc);

select distinct uniprot_id, symbol, residue, site, specified_phosphosite, `avg(logfc_0010)` as logfc, time from table1;

drop table table1_max;
create table table1_max
select distinct uniprot_id, symbol, residue, site, specified_phosphosite, max(`avg(logfc_0010)`) as logfc, group_concat(time separator ";") as 'time'
from table1
group by uniprot_id, symbol, residue, site 
order by abs(max(`avg(logfc_0010)`)) desc;

select distinct * 
from table1_max
where logFC > 1;

drop table table1_min;
create table table1_min
select distinct uniprot_id, symbol, residue, site, specified_phosphosite, min(`avg(logfc_0010)`) as logfc, group_concat(time separator ";") as 'time'
from table1
group by uniprot_id, symbol, residue, site 
order by abs(min(`avg(logfc_0010)`)) desc;

select distinct * 
from table1_min
where logFC < -1;

###################################################################################################
###################################################################################################
########### import and process  top.table.all.csv with specified phospho site for flagging early and late timepoints
########### top.all.sub.csv is the same as top.table.all.csv but with isolated peptide id
### use import wizared in mysql workbench for tables 
### top.all.csv
### specified_phosphosite.csv


### annotate phosphosite
alter table `lipidomics`.`top.all` add column specified_phosphosite int(2);
update `lipidomics`.`top.all` set specified_phosphosite = 0;

update `lipidomics`.`top.all` a
inner join `lipidomics`.`specified_phosphosite` b
on a.ID = b.ID
set a.specified_phosphosite = 1;



### annotate PKA targets selected by mita
alter table `lipidomics`.`top.all` add column interesting_proteins_mita int(2);
update `lipidomics`.`top.all` set interesting_proteins_mita = 0;

update `lipidomics`.`top.all` a
inner join `lipidomics`.`top.pka.mita` b
on a.uniprot_id_10 = b.uniprot_id_10
set a.interesting_proteins_mita = 1;



### annoate description
alter table `lipidomics`.`top.all` add column description mediumtext;

update `lipidomics`.`top.all` a
inner join phosphop b
on a.uniprot_id_10 = b.Uniprot_id
set a.description = b.Description;



### annotate timepoint regulation
alter table `lipidomics`.`top.all` add column early_response int(2);
alter table `lipidomics`.`top.all` add column late_response int(2);

update `lipidomics`.`top.all` set early_response = 0;
update `lipidomics`.`top.all` set late_response = 0;

update `lipidomics`.`top.all` set early_response = 1 
where `adj.p.value_0010` <= 0.05
or `adj.p.value_0030` <= 0.05
or `adj.p.value_0060` <= 0.05
or `adj.p.value_0300` <= 0.05;


update `lipidomics`.`top.all` set late_response = 1 
where `adj.p.value_0600` <= 0.05
or `adj.p.value_0900` <= 0.05
or `adj.p.value_1800` <= 0.05;



### test different hypothesis e.g. time point regulation, with clear site, already specified phosphosite, 
SELECT distinct ID,uniprot_id_10 as uniport_id,description,symbol_0010 as symbol ,Residue,Site,Sequence,specified_phosphosite,
logfc_0010,logfc_0030,logfc_0060,logfc_0300,logfc_0600,logfc_0900,logfc_1800,
`adj.p.value_0010`,`adj.p.value_0030`,`adj.p.value_0060`,`adj.p.value_0300`,`adj.p.value_0600`,`adj.p.value_0900`,`adj.p.value_1800`
FROM lipidomics.`top.all`
where early_response = 1
or late_response = 1)


SELECT distinct ID,uniprot_id_10 as uniport_id,description,symbol_0010 as symbol ,Residue,Site,Sequence,specified_phosphosite,
logfc_0010,logfc_0030,logfc_0060,logfc_0300,logfc_0600,logfc_0900,logfc_1800,
`adj.p.value_0010`,`adj.p.value_0030`,`adj.p.value_0060`,`adj.p.value_0300`,`adj.p.value_0600`,`adj.p.value_0900`,`adj.p.value_1800`
FROM lipidomics.`top.all`
where (early_response = 1
or late_response = 1)
and specified_phosphosite = 1
and Residue != ""

order by abs(logfc_1800) desc;

and abs(logfc_1800) > 0.3
order by `adj.p.value_1800` asc
#and symbol != "RAD23B"


###################################################################################################
###################################################################################################
########### export for heatmap for results

### neutral selection of targets only significant 
SELECT distinct ID,uniprot_id_10 as uniprot_id, description, symbol_0010 as symbol, Residue, Site, Sequence, 
specified_phosphosite, interesting_proteins_mita,
early_response,late_response,
logfc_0010,logfc_0030,logfc_0060,logfc_0300,logfc_0600,logfc_0900,logfc_1800, 
`adj.p.value_0010`,`adj.p.value_0030`,`adj.p.value_0060`,`adj.p.value_0300`,`adj.p.value_0600`,`adj.p.value_0900`,`adj.p.value_1800`
FROM lipidomics.`top.all`
where  (late_response = 1 or late_response = 1) and Residue != ""))
or specified_phosphosite = 1
group by uniprot_id, symbol, residue, site
order by symbol asc;

### select specified phosphosites only but also non significant
SELECT distinct ID,uniprot_id_10 as uniprot_id, description, symbol_0010 as symbol, Residue, Site, Sequence, 
specified_phosphosite, interesting_proteins_mita,
early_response,late_response,
logfc_0010,logfc_0030,logfc_0060,logfc_0300,logfc_0600,logfc_0900,logfc_1800, 
`adj.p.value_0010`,`adj.p.value_0030`,`adj.p.value_0060`,`adj.p.value_0300`,`adj.p.value_0600`,`adj.p.value_0900`,`adj.p.value_1800`
FROM lipidomics.`top.all`
where  specified_phosphosite = 1
group by uniprot_id, symbol, residue, site
order by symbol asc;

drop table heatmap_sepcified_phosphosite;
create table heatmap_sepcified_phosphosite
SELECT distinct a.uniprot_id_10 as uniprot_id, a.symbol_0010 as symbol, a.residue, a.site, a.ID,
avg(b.X0000_1), avg(b.X0000_2), avg(b.X0000_3), avg(b.X0000_4), avg(b.X0000_5), avg(b.X0000_7), avg(b.X0000_8),
avg(b.X0010_1), avg(b.X0010_2), avg(b.X0010_3), avg(b.X0010_5), avg(b.X0010_7), avg(b.X0010_8),
avg(b.X0030_1), avg(b.X0030_2), avg(b.X0030_3), avg(b.X0030_4), avg(b.X0030_5), avg(b.X0030_7), avg(b.X0030_8),
avg(b.X0060_1), avg(b.X0060_2), avg(b.X0060_3), avg(b.X0060_4), avg(b.X0060_5), avg(b.X0060_7), avg(b.X0060_8),
avg(b.X0300_1), avg(b.X0300_2), avg(b.X0300_3), avg(b.X0300_4), avg(b.X0300_5), avg(b.X0300_7), avg(b.X0300_8),
avg(b.X0600_1), avg(b.X0600_2), avg(b.X0600_3), avg(b.X0600_5), avg(b.X0600_7), avg(b.X0600_8),
avg(b.X0900_1), avg(b.X0900_2), avg(b.X0900_3), avg(b.X0900_4), avg(b.X0900_5), avg(b.X0900_7), avg(b.X0900_8),
avg(b.X1800_1), avg(b.X1800_2), avg(b.X1800_3), avg(b.X1800_5), avg(b.X1800_7), avg(b.X1800_8)
FROM lipidomics.`top.all` a
inner join lipidomics.norm_abundance b
on a.ID = b.ID
where  specified_phosphosite = 1
#and (early_response = 1 or middle_response= 1 or late_response = 1)
and residue != ""
group by uniprot_id, symbol, a.residue, a.site
order by symbol asc;

select distinct uniprot_id, symbol, residue, site from heatmap_sepcified_phosphosite;
select distinct uniprot_id, symbol, residue, site from table1;
###################################################################################################
###################################################################################################
########### export all information in one table
########### add petide info to all tables 
########### use the mysql workbecnh import wizard to import the table A08_pp_info

#export for R phosR
use lipidomics; 
select distinct concat(uniprot_id, ";",  gene_symbol, ";", phosphosite, ";", sequence, ";", peptide_id) as peptide_id,
1800_1, 0900_1, 0600_1, 0300_1, 0060_1, 0030_1, 0010_1, 0000_1, 0000_2, 0010_2, 0030_2, 0060_2, 
0300_2, 0600_2, 0900_2, 1800_2, 1800_3, 0600_3, 0060_3, 0010_3, 0900_3, 0300_3, 0030_3, 0000_3, 
0000_4, 0030_4, 0300_4, 0900_4, 0060_4, 0900_5, 0300_5, 0030_5, 0000_5, 1800_5, 0600_5, 0060_5,
0010_5, 0300_7, 0600_7, 0900_7, 1800_7, 0060_7, 0030_7, 0010_7, 0000_7, 0060_8, 0030_8, 0010_8, 
0000_8, 0300_8, 0600_8, 0900_8, 1800_8
  from lipidomics.a08_phosr_v3
  where gene_symbol is not null; 


########### export all phosphopetide with flags and info
SELECT distinct a.ID, a.peptide_id, a.uniprot_id_10 as uniprot_id, b.uniprot_id as uniprot_id_original, a.symbol_0010 as symbol, a.description, a.Residue, a.Site, a.Sequence, 
a.specified_phosphosite, a.interesting_proteins_mita, a.early_response, a.late_response,
b.Modifications, b.`Modification Pattern`, b.isoform, b.`Protein Groups`, b.Proteins, b.PSMs, b.`Positions in Master Proteins`,
b.`Modifications in Master Proteins`,
b.missed_cleavages, b.mhplus, b.`q-value`, b.pep, b.ion_score, b.apex,
logfc_0010,logfc_0030,logfc_0060,logfc_0300,logfc_0600,logfc_0900,logfc_1800,
`adj.p.value_0010`,`adj.p.value_0030`,`adj.p.value_0060`,`adj.p.value_0300`,`adj.p.value_0600`,`adj.p.value_0900`,`adj.p.value_1800`
FROM lipidomics.`top.all` a
left outer join lipidomics.`a08_pp_info` b
on a.peptide_id = b.peptide_id
where a.Residue != "";



###################################################################################################
###################################################################################################
########### top.table for and for GO enrichment analysis
########### collapsed normalised abundance for GSEA analysis

select 
UniprotID as uniprot_id, 
avg(X0000_1), avg(X0000_2), avg(X0000_3), avg(X0000_4), avg(X0000_5), avg(X0000_7), avg(X0000_8),
avg(X0010_1), avg(X0010_2), avg(X0010_3), avg(X0010_5), avg(X0010_7), avg(X0010_8X0010_8), 
avg(X0030_1), avg(X0030_2), avg(X0030_3), avg(X0030_4), avg(X0030_5), avg(X0030_7), avg(X0030_8),
avg(X0060_1), avg(X0060_2), avg(X0060_3), avg(X0060_4), avg(X0060_5), avg(X0060_7), avg(X0060_8),
avg(X0300_1), avg(X0300_2), avg(X0300_3), avg(X0300_4), avg(X0300_5), avg(X0300_7X), avg(X0300_8),
avg(X0600_1), avg(X0600_2), avg(X0600_3), avg(X0600_5), avg(X0600_7), avg(X0600_8),
avg(X0900_1), avg(X0900_2), avg(X0900_3), avg(X0900_4), avg(X0900_5), avg(X0900_7), avg(X0900_8),
avg(X1800_1), avg(X1800_2), avg(X1800_3), avg(X1800_5), avg(X1800_7), avg(X1800_8)
from norm_abundance
group by UniprotID;


