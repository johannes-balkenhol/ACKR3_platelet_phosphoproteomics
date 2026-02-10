# use import wizard for uploading lipidomics.a08_phosr_v3


drop table lipidomics.a08_phosr_v3;
create table lipidomics.a08_phosr_v3
select a.gene_symbol, b.*
from lipidomics.0021_uniprot_mapping a
right join lipidomics.a08_phosr_v2 b
on a.uniprot_id = b.uniprot_id
order by gene_symbol asc;


select count(distinct uniprot_id) from lipidomics.a08_phosr_v3  where gene_symbol is null; 

select distinct uniprot_id from lipidomics.a08_phosr_v3  where gene_symbol is null;

## map gene symbols online in online uniprot_mapping

update lipidomics.a08_phosr_v3 a
inner join lipidomics.mapping_gene_symbol b
on a.uniprot_id = b.uniprot_id
set a.gene_symbol = b.gene_symbol
where a.gene_symbol is null;

select distinct uniprot_id from lipidomics.a08_phosr_v3  where gene_symbol is null;

select distinct uniprot_id from lipidomics.a08_phosr_v3  where gene_symbol is null;

#export for R phosR
select distinct * from lipidomics.a08_phosr_v3; 

select distinct concat(uniprot_id, ";",  gene_symbol, ";", phosphosite, ";", sequence, peptide_id) as peptide_id,
1800_1, 0900_1, 0600_1, 0300_1, 0060_1, 0030_1, 0010_1, 0000_1, 0000_2, 0010_2, 0030_2, 0060_2, 
0300_2, 0600_2, 0900_2, 1800_2, 1800_3, 0600_3, 0060_3, 0010_3, 0900_3, 0300_3, 0030_3, 0000_3, 
0000_4, 0030_4, 0300_4, 0900_4, 0060_4, 0900_5, 0300_5, 0030_5, 0000_5, 1800_5, 0600_5, 0060_5,
0010_5, 0300_7, 0600_7, 0900_7, 1800_7, 0060_7, 0030_7, 0010_7, 0000_7, 0060_8, 0030_8, 0010_8, 
0000_8, 0300_8, 0600_8, 0900_8, 1800_8
  from lipidomics.a08_phosr_v3
  where gene_symbol is not null; 
