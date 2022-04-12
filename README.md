# RBackOffNgram
Back-off n-gram language model for R

This package provides a back-off n-gram for R.

## Example

```{r}
> library(RBackOffNgram)
> x <- "The quick brown fox jumps over the lazy dog"
> words <- strsplit(x," ")[[1]]
> ngram <- ngram_train(words,3,unk_thres=0)
> ngram_lprob(ngram,words)
[1] -7.977968
```
