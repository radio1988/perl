package READ_SYMBOL;

#for *.JGI_symbol file in the PAS-Seq analysis

sub SYMBOL{
        my$line = shift;
	$line =~ s/[\r\n]//;
        my@ele = split /\t/,$line;
	my%pas;
        (
	$pas{scaffold},
	$pas{direction},
	$pas{position},
	$pas{count},
	$pas{read_id},
	$pas{cigar},
	$pas{cigar_pos},
	$pas{mq},
	$pas{read_seq},
	$pas{seq100},
	$pas{gene_acc},#gene acc
	$pas{class_code},
	$pas{read_id},#read id
	$pas{read_id2},#read id
	#$pas{gene_acc_bug},#gene acc when mapped, when not mapped, gives -
	#$pas{gene_page},
	#$pas{gene_symbol},
	#$pas{gene_desc}
	
	
        )=@ele;
	#print"line:$line;\n";
	#print"symbol:$pas{gene_symbol}, desc$pas{gene_desc};\n";
	#sleep 2;
        return %pas;


}
1;

