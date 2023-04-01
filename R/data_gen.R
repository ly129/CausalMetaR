generate_simulated_data <- function(n = 2000) {
  prevmi <- rbinom(n, size = 1, prob = 0.50)
  age <- runif(n)
  angina <- rbinom(n, size = 1, prob = 0.50)
  ejecfr <- runif(n)
  lvscor <- runif(n)
  prxlad31 <- runif(n)
  prxves31 <- runif(n)

  treatment <- rbinom(n, 1, 0.5)

  source1 <- 1
  source2 <- prevmi + age + angina + ejecfr
  source3 <- -1 + 0.5*(prevmi + age + angina + ejecfr)
  vProb <- exp(cbind(source1, source2, source3))
  mChoices <- apply(vProb, 1, rmultinom, n = 1, size = 1)

  outcome <- rbinom(n, 1, plogis(-1 + 0.5*(prevmi + age + angina + ejecfr) * treatment))

  return(data.frame(prevmi,
                    age,
                    angina,
                    ejecfr,
                    lvscor,
                    prxlad31,
                    prxves31,
                    treatment,
                    source1 = mChoices[1,],
                    source2 = mChoices[2,],
                    source3 = mChoices[3,],
                    outcome))
}
