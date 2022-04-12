#library(hash)
#library(R6)

Ngram_node <- R6::R6Class("Ngram_node",
  public=list(
    rvocab=NULL,
    vocab=NULL,
    subnode=NULL,
    count=NULL,
    total=0,
    alpha=NULL,
    beta=NULL,
    training=TRUE,
    initialize=function() {
      self$rvocab <- hash::hash("<UNK>",1)
      self$vocab <- c("<UNK>")
      self$subnode <- list(FALSE)
      self$count <- c(0)
      self$total <- 0
      self$alpha <- NULL
      self$beta <- NULL
    },
    put=function(words,cnt=1) {
      w <- words[1]
      if (!hash::has.key(w,self$rvocab)) {
        wpos <- length(self$rvocab)+1
        self$rvocab[[w]] <- wpos
        self$vocab[wpos] <- w
        self$count[wpos] <- 0
        self$subnode[[wpos]] <- FALSE
      } else {
        wpos <- self$rvocab[[w]]
      }
      self$count[wpos] <- self$count[wpos]+cnt
      self$total <- self$total+cnt
      if (length(words) > 1) {
        nextnode <- self$subnode[[wpos]]
        if (is.logical(nextnode)) {
          nextnode <- Ngram_node$new()
          self$subnode[[wpos]] <- nextnode
        }
        nextnode$put(words[2:length(words)],cnt)
      }
    },
    getCount=function(words) {
      words <- self$wordReplace(words)
      self$getCount0(words)
    },
    getCount0=function(words) {
      L <- length(words)
      w <- words[1]
      if (hash::has.key(w,self$rvocab)) {
        wpos <- self$rvocab[[w]]
      } else {
        return(0)
      }
      if (L == 1) {
        return(self$count[wpos])
      }
      sbn <- self$subnode[[wpos]]
      if (is.logical(sbn)) return(0)
      sbn$getCount0(words[2:L])
    },
    wordReplace=function(words,func) {
      ws <- sapply(words,function(w) {
        if (hash::has.key(w,self$rvocab)) {
          return(w)
        }
        "<UNK>"
      })
      if (!self$training && !self$has_unk() && any(ws %in% "<UNK>")) {
        print(ws)
        stop("Unknown words are not accepted for models without UNK")
      }
      ws
    },
    getNext=function(words) {
      words <- self$wordReplace(words)
      self$getNext0(words)
    },
    getNext0=function(words) {
      if (length(words) == 0) {
        return(data.frame(word=self$vocab,count=self$count))
      }
      wpos <- self$rvocab[[words[1]]]
      if (length(words) > 1) {
        return(self$subnode[[wpos]]$getNext0(words[2:length(words)]))
      }
      self$subnode[[wpos]]$getNext0(c())
    },
    mlprob=function(words) {
      words <- self$wordReplace(words)
      L <- length(words)
      if (L == 1) {
        wpos <- self$rvocab[[words]]
        return(self$count[wpos]/self$total)
      }
      n1 <- self$getCount(words)
      n2 <- self$getCount(words[1:(length(words)-1)])
      n1/n2
    },
    calc_discount=function() {
      self$calc_discount0()
      self$calc_discount1(self,NULL)
    },
    calc_discount0=function() {
      self$alpha <- rep(1,length(self$vocab))
      self$beta <- rep(1,length(self$vocab))
      for (i in 1:length(self$vocab)) {
        if (is.logical(self$subnode[[i]])) {
          self$alpha[i] <- 0
        } else {
          r <- length(self$subnode[[i]]$vocab)
          n <- self$count[i]
          self$alpha[i] <- 1-r/(n+r)
          self$subnode[[i]]$calc_discount0()
        }
      }
    },
    calc_discount1=function(top,words) {
      for (i in 1:length(self$vocab)) {
        total <- 0
        ww <- c(words,self$vocab[i])
        L <- length(ww)
        sbn <- self$subnode[[i]]
        if (is.logical(sbn)) next
        for (j in 1:length(sbn$vocab)) {
          if (L == 1) {
            # unigram
            total <- total+top$boprob(sbn$vocab[j])
          } else {
            # bigram or more
            total <- total+top$boprob(c(ww[(L+1):length(ww)],sbn$vocab[j]))
          }
        }
        self$beta[i] <- 1/(1-total)
        sbn$calc_discount1(top,ww)
      }
    },
    getAlpha=function(words) {
      words <- self$wordReplace(words)
      self$getAlpha0(words)
    },
    getAlpha0=function(words) {
      L <- length(words)
      w <- words[1]
      wpos <- self$rvocab[[w]]
      if (L == 1) {
        return(c(self$alpha[wpos],self$beta[wpos]))
      }
      sbn <- self$subnode[[wpos]]
      if (is.logical(sbn)) return(c(1,1))
      sbn$getAlpha0(words[2:L])
    },
    boprob=function(words) {
      words <- self$wordReplace(words)
      self$boprob0(words)
    },
    boprob0=function(words) {
      L <- length(words)
      if (L == 1) {
        return(self$getCount(words)/self$total)
      }
      n1 <- self$getCount(words)
      n2 <- self$getCount(words[1:(length(words)-1)])
      a <- self$getAlpha0(words[1:(length(words)-1)])
      if (n1 > 0) {
        return(a[1]*n1/n2)
      }
      if (n2 > 0) {
        return((1-a[1])*a[2]*self$boprob0(words[2:length(words)]))
      }
      self$boprob0(words[2:length(words)])
    },
    has_unk=function() {
      self$count[1] > 0
    }
  )
)

#' Train back-off n-gram language model
#' @param x a vector of characters (words)
#' @param n length of n-gram
#' @param unk_thres threshold of word frequency.
#'    If the frequency of the word is no more than this value, that word is replaced with UNK symbol.
#' @param start_sym if not NA, this symbol is appended at the top of x
#' @param end_sym if not NA, this symbol is appended at the end of x
#' @returns n-gram model
#' @importFrom hash hash has.key
#' @export
ngram_train <- function(x,n,unk_thres=1,start_sym=NA,end_sym=NA) {
  wfreq <- hash::hash()
  for (i in 1:length(x)) {
    if (hash::has.key(x[i],wfreq)) {
      wfreq[[x[i]]] <- wfreq[[x[i]]]+1
    } else {
      wfreq[[x[i]]] <- 1
    }
  }
  for (i in 1:length(x)) {
    if (wfreq[[x[i]]] <= unk_thres) {
      x[i] <- "<UNK>"
    }
  }
  if (!is.na(start_sym)) {
    x <- c(rep(start_sym,n-1),x)
  }
  if (!is.na(end_sym)) {
    x <- c(x,end_sym)
  }
  node <- Ngram_node$new()
  if (n > 1) {
    L <- length(x)
    for (i in 1:(n-1)) {
      node$put(x[1:i])
      node$put(x[(L-i+1):L])
    }
  }
  for (i in n:length(x)) {
    node$put(x[(i-n+1):i])
  }
  node$calc_discount()
  node$training <- FALSE
  list(n=n,node=node)
}

#' Calculate log probability of the sequence
#' @param ngram An n-gram model
#' @param x q sequence of characters (words). It should be longer than the length of the n-gram.
#' @return log probability of x
#' @export
ngram_lprob <- function(ngram,x) {
  total_lprob <- 0
  n <- ngram$n
  if (length(x) < n) {
    stop("Input length should be longer than n-gram length")
  }
  for (i in n:length(x)) {
    total_lprob <- total_lprob+
      log(ngram$node$boprob(x[(i-n+1):i]))
  }
  total_lprob
}


