#A few examples of get_blast_seeds

get_blast_seeds("TAGAACAGGCTCCTCTAG" , "TTAGATACCCCACTATGC" ,
                "D:/blast_seeds", "12S_V5F1", "D:/accessionTaxa.sql",
                organism = c("7776"), MAX_TARGET_PER_TEMPLATE = 10)
