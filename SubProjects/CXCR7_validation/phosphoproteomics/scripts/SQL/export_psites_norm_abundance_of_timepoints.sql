#################################################################
##### chosse database
use lipidomics; 


#################################################################
##### select significant norm abundance of different timepoints

## early timepoint
select distinct a.*
from norm_abundance a
inner join 
(select distinct ID from a08_all_response_all_targets_all_20122021 
where `adj.p.value_0010` < 0.06 or `adj.p.value_0030` < 0.06 or `adj.p.value_0060` < 0.06  ) b
on a.ID = b.ID; 

## intermediate timepoint
select distinct a.*
from norm_abundance a
inner join 
(select distinct ID from a08_all_response_all_targets_all_20122021 
where `intermediate responder 300+600 sec` = 1) b
on a.ID = b.ID; 

## late timepoint
select distinct a.*
from norm_abundance a
inner join 
(select distinct ID from a08_all_response_all_targets_all_20122021 
where `late responder 900 + 1800 sec` = 1) b
on a.ID = b.ID;


#################################################################
#### select significant log2fc of different timepoints

## early timepoint
select distinct a.ID, logfc_0010, logfc_0030, logfc_0060, logfc_0300, logfc_0600, logfc_0900, logfc_1800
from norm_abundance a
inner join 
(select distinct * from a08_all_response_all_targets_all_20122021 
where `adj.p.value_0010` < 0.06 or `adj.p.value_0030` < 0.06 or `adj.p.value_0060` < 0.06  ) b
on a.ID = b.ID; 

## intermediate timepoint
select distinct a.ID, logfc_0010, logfc_0030, logfc_0060, logfc_0300, logfc_0600, logfc_0900, logfc_1800
from norm_abundance a
inner join 
(select distinct * from a08_all_response_all_targets_all_20122021 
where `intermediate responder 300+600 sec` = 1) b
on a.ID = b.ID; 

## late timepoint
select distinct a.ID, logfc_0010, logfc_0030, logfc_0060, logfc_0300, logfc_0600, logfc_0900, logfc_1800
from norm_abundance a
inner join 
(select distinct * from a08_all_response_all_targets_all_20122021 
where `late responder 900 + 1800 sec` = 1) b
on a.ID = b.ID;

#################################################################
##### averaging psite
select distinct a.ID, 
avg(X0000_1), avg(X0000_2), avg(X0000_3), avg(X0000_4), avg(X0000_5), avg(X0000_7), avg(X0000_8), 
avg(X0010_1), avg(X0010_2), avg(X0010_3), avg(X0010_5), avg(X0010_7), avg(X0010_8), 
avg(X0030_1), avg(X0030_2), avg(X0030_3), avg(X0030_4), avg(X0030_5), avg(X0030_7), avg(X0030_8), 
avg(X0060_1), avg(X0060_2), avg(X0060_3), avg(X0060_4), avg(X0060_5), avg(X0060_7), avg(X0060_8), 
avg(X0300_1), avg(X0300_2), avg(X0300_3), avg(X0300_4), avg(X0300_5), avg(X0300_7), avg(X0300_8), 
avg(X0600_1), avg(X0600_2), avg(X0600_3), avg(X0600_5), avg(X0600_7), avg(X0600_8), 
avg(X0900_1), avg(X0900_2), avg(X0900_3), avg(X0900_4), avg(X0900_5), avg(X0900_7), avg(X0900_8), 
avg(X1800_1), avg(X1800_2), avg(X1800_3), avg(X1800_5), avg(X1800_7), avg(X1800_8)
from norm_abundance a
inner join 
(select distinct * from a08_all_response_all_targets_all_20122021 
where  `immediate early 10 sec` = 1  or
`early responder 30  + 60 sec` = 1 
)
 b
on a.ID = b.ID
group by symbol, Site, Residue; 


#################################################################
#### filter here duplicate psites
## late timepoint
select distinct a.*
from norm_abundance a
inner join 
(select distinct * from a08_all_response_all_targets_all_20122021 
where `late responder 900 + 1800 sec` = 1) b
on a.ID = b.ID
group by symbol, Residue, Site;


#################################################################
#### here special kinase regulation

## early timepoint
select distinct a.ID, logfc_0010, logfc_0030, logfc_0060, logfc_0300, logfc_0600, logfc_0900, logfc_1800
from norm_abundance a
inner join 
(select distinct * from a08_all_response_all_targets_all_20122021 
where `adj.p.value_0010` < 0.06 or `adj.p.value_0030` < 0.06 or `adj.p.value_0060` < 0.06  ) b
on a.ID = b.ID
where symbol like "lyst%";
#PTPN12 Interacts with TGFB1I1 
#

## intermediate timepoint
select distinct a.ID, 
logfc_0010, logfc_0030, logfc_0060, logfc_0300, logfc_0600, logfc_0900, logfc_1800,
(X0000_1 + X0000_2 + X0000_3 + X0000_4 + X0000_5 + X0000_7) / 6 as X0000,
(X0010_1 + X0010_2 + X0010_3 + X0010_5 + X0010_7 + X0010_8) / 6 as X0010,
(X0030_1 + X0030_2 + X0030_3 + X0030_4 + X0030_5 + X0030_7 + X0030_8) / 7 as X0030,
(X0060_1 + X0060_2 + X0060_3 + X0060_4 + X0060_5 + X0060_7 + X0060_8) / 7 as X0060,
(X0300_1 + X0300_2 + X0300_3 + X0300_4 + X0300_5 + X0300_7 + X0300_8) / 7 as X0300,
(X0600_1 + X0600_2 + X0600_3 + X0600_5 + X0600_7 + X0600_8) / 6 as X0600,
(X0900_1 + X0900_2 + X0900_3 + X0900_4 + X0900_5 + X0900_7 + X0900_8) / 7 as X9000, 
(X1800_1 + X1800_2 + X1800_3 + X1800_5 + X1800_7 + X1800_8) / 6 as X1800
from raw_abundance a
inner join 
(select distinct * from a08_all_response_all_targets_all_20122021
where `adj.p.value_0300` < 0.05 or `adj.p.value_0600` < 0.05
or `adj.p.value_0010` < 0.06 or `adj.p.value_0030` < 0.06 or `adj.p.value_0060` < 0.06 ) b
on a.ID = b.ID
where symbol like "prka%";
#stk25 stk39
#https://www.uniprot.org/uniprot/Q9UEW8
#https://www.ebi.ac.uk/QuickGO/term/GO:0038146

#PTPN12 Interacts with TGFB1I1 


## late timepoint
select distinct a.ID, logfc_0010, logfc_0030, logfc_0060, logfc_0300, logfc_0600, logfc_0900, logfc_1800
from norm_abundance a
inner join 
(select distinct * from a08_all_response_all_targets_all_20122021 
where `adj.p.value_0900` < 0.08 or `adj.p.value_1800` < 0.08) b
on a.ID = b.ID
where symbol like "PRK%";



#################################################################
## filter early timepoint regulated sites for the string network
select distinct b.uniprot_id
#, a.ID, logfc_0010, logfc_0030, logfc_0060, logfc_0300, logfc_0600, logfc_0900, logfc_1800
from norm_abundance a
inner join 
(select distinct * from a08_all_response_all_targets_all_20122021 
where `adj.p.value_0300` < 0.05 or `adj.p.value_0600` < 0.05 ) b
#where `adj.p.value_0010` < 0.06 or `adj.p.value_0030` < 0.06  or `adj.p.value_0060` < 0.06) b
on a.ID = b.ID;
#PTPN12 Interacts with TGFB1I1 
#CXCR7 interacts with EGFR interacts with PTPN12