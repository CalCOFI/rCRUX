#' Query a pair of primers using NCBI's Primer-BLAST, if primers contain iupac
#'
#' @details 
#' PrimerTree
#' Modified 2022-08-10 by Luna Gal
#' Removed num_aligns argument from primer_search
#' Instead num_targets_with_primers is determined by the user
#' Changed the calculation for duration of blast submission
#'
#' ambiguity codes, enumerate all possible combinations and combine the
#' results.
#'
#' @param forward forward primer to search by 5'-3' on plus strand
#' @param reverse reverse primer to search by 5'-3' on minus strand
#' @param num_permutations the number of primer permutations to search, if the degenerate bases
#'        cause more than this number of permutations to exist, this number will be
#'        sampled from all possible permutations.
#' @param ... additional arguments passed to Primer-Blast, BLAST_primer
#' @param .parallel if 'TRUE', perform in parallel, using parallel backend
#'        provided by foreach
#' @param .progress name of the progress bar to use, see 'create_progress_bar'
#' 
#' @return a list of results for each primer?
primer_search <- function(forward, reverse, num_permutations = 25, ..., .parallel = FALSE, .progress = 'none'){
  
  if (missing(forward) || missing(reverse)) {
    BLAST_primer()
    primers = enumerate_primers(forward, reverse)
    num_primers = nrow(primers)
  }
  
  if (num_primers > num_permutations) {
    warning(immediate.=T, 'primers have ', num_primers, ' possible combinations due to degenerate bases, sampling ', num_permutations, " primers, use 'num_permutations' to change")
    primers = primers[ sample.int(num_primers, num_permutations, replace=F), ]
  }
  
  message('BLASTing ', nrow(primers), ' primer combinations')
  #enumerate all combinations to handle ambiguity codes
  plyr::alply(.data = primers, 
              .margins = 1, 
              .expand = F,
              .parallel = .parallel, 
              .progress = .progress,
              function(row) {
                BLAST_primer(row$forward, row$reverse, ...)
              })
}

enumerate_primers = function(forward, reverse){
  forward_primers = enumerate_ambiguity(forward)
  data.frame(forward=forward_primers,
             reverse=rep(enumerate_ambiguity(reverse),
                         each=length(forward_primers)),
             stringsAsFactors = FALSE)
}

enumerate_ambiguity = function(sequence){
  
  iupac = list( "M" = list("A", "C"),
                "R" = list("A", "G"),
                "W" = list("A", "T"),
                "S" = list("C", "G"),
                "Y" = list("C", "T"),
                "K" = list("G", "T"),
                "V" = list("A", "C", "G"),
                "H" = list("A", "C", "T"),
                "D" = list("A", "G", "T"),
                "B" = list("C", "G", "T"),
                "N" = list("A", "C", "G", "T"),
                "I" = list("A", "T", "C"))
  
  search_regex = paste(names(iupac), collapse='|')
  locs = stringr::str_locate_all(sequence, search_regex)
  sequences = list()
  count = 1
  for (i in seq_len(nrow(locs[[1]]))){
    loc = locs[[1]][i,]
    ambiguity = stringr::str_sub(sequence, loc[1], loc[2])
    for(type in iupac[[ambiguity]]){
      new_seq = sequence
      stringr::str_sub(new_seq, loc[1], loc[2]) <- type
      sequences[[count]] = enumerate_ambiguity(new_seq)
      count = count + 1
    }
    return(unlist(sequences))
  }
  return(sequence)
}

print_options = function(options){
  output = utils::capture.output(print(
    options[ is.na(options$type) | options$type != 'hidden',
             c('name', 'type', 'defval') ])
  )
  
  message(paste(output, "\n", sep=""))
}

#' Query a pair of primers using ncbi's Primer-BLAST
#' 
#' Query a pair of primers using ncbi's Primer-BLAST at https://www.ncbi.nlm.nih.gov/tools/primer-blast/
#' 
#' @param forward forward sequence
#' @param reverse reverse sequence
#' @param ... passed onto the POST request (see NCBI for details)
#' @param organism organism
#' @param primer_specificity_database primer_specificity_database='nt'
#' @param exclude_env exclude_env='on'
#' 
#' @return httr response object of the query, pass to parse_primer_hits to parse the results.
BLAST_primer = function(forward, reverse, ..., organism='',
                        primer_specificity_database='nt', exclude_env='on'){
  
  # modify a function to check the status and retry until success
  retry = function(fun, num_retry=5, ...){
    function(...){
      res = fun(...)
      itr = 0
      status = httr::http_status(res)
      while(itr < num_retry && inherits(res, 'response') && tolower(status$category) != 'success'){
        warning('request failed, retry attempt ', itr+1)
        res = fun(...)
        status = httr::http_status(res)
        itr = itr + 1
        #sleep to avoid hitting NCBI query rate limit :'-(
        Sys.sleep(0.4)
      }
      res
    }
  }
  
  GET_retry = retry(httr::GET)
  POST_retry = retry(httr::POST)
  
  url = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/'
  form = GET_retry(url)
  
  content = parsable_html(form)
  
  all_options = get_options(content)
  
  if(missing(forward) || missing(reverse)){
    print_options(all_options)
    stop('No primers specified')
  }
  
  options = list(..., primer_left_input=forward, primer_right_input=reverse,
                 organism=organism,
                 primer_specificity_database=primer_specificity_database,
                 exclude_env=exclude_env,
                 search_specific_primer='on')
  
  names(options) = toupper(names(options))
  
  match_args = pmatch(names(options), all_options$name)
  bad_args = is.na(match_args)
  
  if(any(bad_args)){
    print_options(all_options)
    stop(paste(names(options)[bad_args], collapse=','), ' not valid option\n')
  }
  
  options = get_defaults(options, all_options)
  
  start_time = lubridate::now()
  
  message('Submitting Primer-BLAST query')
  response = POST_retry(paste(url, 'primertool.cgi', sep=''), body=options)
  
  values = get_refresh_from_meta(response)
  
  while(length(values) > 0){
    message('BLAST alignment processing, refreshing in ', values[1], ' seconds...')
    Sys.sleep(values[1])
    response = GET_retry(values[2])
    
    values = get_refresh_from_meta(response)
  }
  
  end_time <- lubridate::now()
  total_time = end_time - start_time
  message('BLAST alignment completed in ', total_time , ' seconds')
  response
}

#' Parse the primer hits
#'
#' @param response a httr response object obtained from \code{\link{primer_search}}
#' @noRd
parse_primer_hits = function(response){
  content = parsable_html(response)
  plyr::rbind.fill(XML::xpathApply(content, '//pre', parse_pre))
}

parse_a = function(a){
  #links like entrez/viewer.fcgi?db=nucleotide&id=452085006
  m = regexpr('id=\\d+', XML::xmlAttrs(a)['href'])
  gi = gsub('id=', '', unlist(regmatches(XML::xmlAttrs(a)['href'], m)))
  if(length(gi) == 0L){
    #links like nucleotide/449036831?from=1107741&to=1107929&report=gbwithparts
    m = regexpr('nucleotide/\\d+', XML::xmlAttrs(a)['href'])
    gi = gsub('nucleotide/', '', unlist(regmatches(XML::xmlAttrs(a)['href'], m)))
  }
  if(length(gi) == 0L) gi = NA
  data.frame(gi=as.character(gi), accession=as.character(XML::xmlValue(a)))
}

parse_pre = function(pre){
  pre_text = XML::xmlValue(pre)
  
  a = XML::getNodeSet(pre, './preceding-sibling::a[1]')
  if(length(a) <= 0)
    stop('Parsing failed for ', pre_text)
  
  ids = parse_a(a[[1]])
  
  product_length_regex = 'product length = (\\d+)'
  template_regex = 'Template[^\\d]+(\\d+)[^.ACGT]+([.ACGT]+)[^\\d]+(\\d+)'
  full_regex = paste('[\\S\\W]*', product_length_regex, '[\\S\\W]*?',
                     template_regex, '[\\S\\W]*', template_regex, '[\\S\\W]*', sep='')
  values = stringr::str_split(gsub(full_regex, paste('\\', 1:8, sep='', collapse='|'),
                                   pre_text, perl=T), '[|]')[[1]]
  data.frame(
    ids,
    product_length = as.numeric(values[1]),
    mismatch_forward = stringr::str_count(values[3], '[ACGT]'),
    mismatch_reverse = stringr::str_count(values[6], '[ACGT]'),
    forward_start = as.numeric(values[2]),
    forward_stop = as.numeric(values[4]),
    reverse_start = as.numeric(values[5]),
    reverse_stop = as.numeric(values[7]),
    product_start = min(as.numeric(values[c(2, 4, 5, 7)])),
    product_stop = max(as.numeric(values[c(2, 4, 5, 7)]))
  )
}

get_refresh_from_meta = function(response){
  content = parsable_html(response)
  meta = content['//meta[@http-equiv="Refresh"]']
  if(length(meta) > 0){
    values = stringr::str_split(XML::xmlAttrs(meta[[1]])['content'], '; URL=')[[1]]
    return(values)
  }
  return()
}

get_defaults = function(set_options, options){
  #only look at values with set defaults
  options = options[ !is.na(options$defval), ]
  unchanged_options = setdiff(options$name, names(set_options))
  default_values = as.character(options[ match(unchanged_options, options$name), 'defval' ])
  names(default_values) = unchanged_options
  c(set_options, default_values)
}

get_options = function(content){
  options = plyr::rbind.fill(XML::xpathApply(content, '//form//input | //form//select', parse_attributes))
  options$type = as.character(options$type)
  
  #add dropdown type if they are NA
  options$type[is.na(options$type)] <- 'dropdown'
  
  #make default values for checkboxes on or off rather than checked or unchecked
  options$defval = as.character(options$defval)
  check_map = c('checked' = 'on', 'unchecked' = '')
  checkboxes = which(options$type == 'checkbox')
  options$defval[ checkboxes ] = check_map[ options$defval[ checkboxes ] ]
  
  options[ options$type != 'hidden',  c('name', 'type', 'defval') ]
}

parse_attributes = function(x){
  as.data.frame(t(XML::xmlAttrs(x)))
}

parsable_html = function(response){
  txt <- httr::content(response, as='text', encoding='UTF-8')
  
  # remove any unicode characters
  Encoding(txt) <- "UTF-8"
  txt <- iconv(txt, "UTF-8", "ASCII", sub = "")
  
  #this gsub regex is to remove the definition lines, some of which have
  #  bracketed <junk> in them, which messes up the parsing
  txt <- gsub('("new_entrez".*?</a>).*?<pre>\n\n', '\\1\n<pre>', txt)
  XML::htmlParse(txt)
}

filter_duplicates = function(hits){
  location_columns = c('accession', 'forward_start', 'forward_stop', 'reverse_start', 'reverse_stop')
  hits[!duplicated(t(apply(hits[location_columns], 1, range))),]
}
