#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

###############	open directory and extratct files ending in .graphml
my $dir = './';
my @files;
my %protein_list_complete;
my %edge_list_complete;
my %edges_help;
my %edge_list_h;
my %species_list_h;
my $node_id_global=0;

my @reactions_sort;
my @protlist_sort;
my $edgecount=0;
my $node_id;
my $nodecount_firstSet=0;


opendir(DIR, $dir) or die $!;
while (my $file = readdir(DIR)) {
  # Use a regular expression to find files ending in .graphml
  if ($file =~ m/\.graphml$/){
    push (@files, $file);
  }
}
closedir(DIR);


#############################################################################
##	building global protein List
#############################################################################
foreach my $f(@files){
open (INPUT_FILE, "./$f") or die $1;
my %protein_list_local;
my $node_id_local;

  while (my $graphml=<INPUT_FILE>){

    if ($graphml=~/<node\sid\=\"n(\d+)\">/){
      $node_id_local=$1;
	  
    }

    #elsif ($graphml=~/<y:NodeLabel.+?>(.+)?<y:LabelModel>$/){
	elsif ($graphml=~/<y:NodeLabel.+?>(.+)<y:LabelModel>?/){
	#print "hallo".$1."\n";
      $protein_list_local{$node_id_local}=$1;  
	        if (defined $protein_list_complete{$1}){}
      else{
	$node_id_global++;
	$protein_list_complete{$1}=$node_id_global;
      }
    }
  }
close (INPUT_FILE) or die $!;
}


#############################################################################
##	building global edge List
#############################################################################
foreach my $f(@files){
	open (INPUT_FILE, "./$f") or die $1;
	my %protein_list_local;
	my $node_id_local;
	my $edge_name;
	my @reac_type_a;

	#print Dumper(%protein_list_complete);

	while (my $graphml=<INPUT_FILE>){
		if ($graphml=~/<node\sid\=\"n(\d+)\">/){
		  $node_id_local=$1;
		}
		elsif ($graphml=~/<y:NodeLabel.+?>(.+)?<y:LabelModel>?/){
		  $protein_list_local{$node_id_local}=$1;
		  #print "NodeID   ".$node_id_local."    NodeName   ".$1."\n";
		}
	}
	close (INPUT_FILE) or die $!;
	
	print Dumper(\%protein_list_local);
	
	open (INPUT_FILE, "./$f") or die $1;
	while (my $graphml=<INPUT_FILE>){
		if ($graphml=~/<edge\sid="e\d+"\ssource="n(\d+)"\starget="n(\d+)">/){
		  print "NodeID1   ".$1."    NodeID2   ".$2."\n";
		  print "NodeName1   ".$protein_list_local{$1}."    NodeName2   ".$protein_list_local{$2}."\n";
		  $edge_name=$protein_list_local{$1}."_".$protein_list_local{$2};
			   if (defined $edge_list_complete{$edge_name}){}
		  else {
		@{$edge_list_complete{$edge_name}}=($protein_list_complete{$protein_list_local{$1}},$protein_list_complete{$protein_list_local{$2}});
		  }
		}
		elsif ($graphml=~/<y:Arrows\ssource="\w+?"\starget="(\w+)"\/>$/) {
		  if ($1 eq "standard"){
		push( @{ $edge_list_complete{$edge_name} }, ("standard","#008000")); 
		  }
		  else{
		push( @{ $edge_list_complete{$edge_name} },("t_shape","#FF0000"));
		  }
		}
	}
	
	#print Dumper(\%edge_list_local);
	
	close (INPUT_FILE) or die $!;
}


@reactions_sort = (keys %edge_list_complete);

@protlist_sort = sort { $protein_list_complete{$a} <=>  $protein_list_complete{$b} } (keys %protein_list_complete);

#print Dumper(\%protein_list_complete);
################# for TXT output ###############
my %rev_protlist;
foreach my $key (keys %protein_list_complete){
  $rev_protlist{$protein_list_complete{$key}} = $key;
}	
my %edgetypes;
$edgetypes{"standard"} = "activation";
$edgetypes{"t_shape"} = "inhibition";
################################################
# print Dumper(\%protein_list_complete);
# print Dumper(\%rev_protlist);
# print Dumper(\%edge_list_complete);
open (TXT, ">./combined.txt") or die $!;
# print TXT Dumper(\%edge_list_complete);
Reac_sub_txt();
close TXT;
open (OUT, ">./combined.graphml") or die $!;


header();
Prot_sub();
Reac_sub();
close_out();
close OUT or die $!;
foreach my $fi (@files){
  print "\"$fi\"\t ";
}
print "\nwhere integrated into \"../combined.graphml\".\n\n";


###############################
##	subroutines for 2graphml_to_1graphml-converter
###############################
sub header {
print OUT '<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:y="http://www.yworks.com/xml/graphml" xmlns:yed="http://www.yworks.com/xml/yed/3" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd">
  <!--Created by yFiles for Java 2.11-->
  <key for="graphml" id="d0" yfiles.type="resources"/>
  <key for="port" id="d1" yfiles.type="portgraphics"/>
  <key for="port" id="d2" yfiles.type="portgeometry"/>
  <key for="port" id="d3" yfiles.type="portuserdata"/>
  <key attr.name="url" attr.type="string" for="node" id="d4"/>
  <key attr.name="description" attr.type="string" for="node" id="d5"/>
  <key for="node" id="d6" yfiles.type="nodegraphics"/>
  <key attr.name="Beschreibung" attr.type="string" for="graph" id="d7"/>
  <key attr.name="url" attr.type="string" for="edge" id="d8"/>
  <key attr.name="description" attr.type="string" for="edge" id="d9"/>
  <key for="edge" id="d10" yfiles.type="edgegraphics"/>
  <graph edgedefault="directed" id="G">
    <data key="d7"/>';
}

sub Prot_sub {  

foreach my $prot_key(@protlist_sort)
{
print OUT '<node id="';
######################################

print OUT 'n'.$protein_list_complete{$prot_key}.'">';


  #######################################

print OUT '<data key="d6">
        <y:ShapeNode>
          <y:Geometry height="45.0" width="45.0" x="275.0" y="271.0"/>
          <y:Fill hasColor="false" transparent="false"/>
          <y:BorderStyle color="#000000" type="line" width="1.0"/>
          <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="16" fontStyle="bold" hasBackgroundColor="false" hasLineColor="false" height="17.96875" modelName="custom" textColor="#000000" visible="true" width="26.88671875" x="1.556640625" y="6.015625">'.$prot_key.'<y:LabelModel>
              <y:SmartNodeLabelModel distance="4.0"/>
            </y:LabelModel>
            <y:ModelParameter>
              <y:SmartNodeLabelModelParameter labelRatioX="0.0" labelRatioY="0.0" nodeRatioX="0.0" nodeRatioY="0.0" offsetX="0.0" offsetY="0.0" upX="0.0" upY="-1.0"/>
            </y:ModelParameter>
          </y:NodeLabel>
          <y:Shape type="ellipse"/>
        </y:ShapeNode>
      </data>
    </node>';
  }
}

sub Reac_sub {
foreach my $key_r(@reactions_sort)
  {
  $edgecount++;
  my @re_data=@{$edge_list_complete{$key_r}};



  if (defined $re_data[2])
    {

    print OUT '<edge id="e'.$edgecount.'" source="n'.$re_data[0].'" target="n'.$re_data[1].'">
      <data key="d9"/>
      <data key="d10">
        <y:PolyLineEdge>
          <y:Path sx="0.0" sy="0.0" tx="0.0" ty="0.0">
            <y:Point x="552.0" y="365.0"/>
          </y:Path>
          <y:LineStyle color="'.$re_data[3].'" type="line" width="1.0"/>
          <y:Arrows source="none" target="'.$re_data[2].'"/>
          <y:BendStyle smoothed="false"/>
        </y:PolyLineEdge>
      </data>
    </edge>';
      }
    }
}

sub close_out {
print OUT '  </graph>
  <data key="d0">
    <y:Resources/>
  </data>
</graphml>';

}

############### for TXT output ############
sub Reac_sub_txt {
foreach my $key_r(@reactions_sort)
  {
  $edgecount++;
  my @re_data=@{$edge_list_complete{$key_r}};



  if (defined $re_data[2])
    {
    print TXT $rev_protlist{$re_data[0]}."\t".$edgetypes{$re_data[2]}."\t".$rev_protlist{$re_data[1]}."\n";

      }
    }
}
