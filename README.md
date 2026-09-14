# Asset Bubble Detection in the Magnificent Seven

## An Econometric Approach Using PSY/GSADF Bubble Detection and Near-Unit-Root Inference

This repository contains the data and R implementation developed for my bachelor's thesis at the University of Trieste.

The project investigates whether the rapid appreciation of the Magnificent Seven technology stocks reflects fundamental growth, broader market movements, or idiosyncratic speculative exuberance.

The empirical framework combines recursive right-tailed unit-root testing with increasingly strict controls for fundamentals and systematic risk.

---

## Research Question

Do the price dynamics of the Magnificent Seven exhibit statistically detectable speculative bubbles after controlling for firm fundamentals, broad-market movements, and systematic equity risk factors?

---

## Assets Analyzed

- Apple (AAPL)
- Microsoft (MSFT)
- Alphabet (GOOG/GOOGL)
- Amazon (AMZN)
- NVIDIA (NVDA)
- Meta Platforms (META)
- Tesla (TSLA)

---

## Empirical Specifications

Four alternative specifications are used.

### 1. Price-to-Sales Ratio

Controls for firm-level fundamental growth by comparing stock prices with company revenues.

### 2. Price-to-S&P 500 Ratio

Controls for broad U.S. equity-market movements and tests whether individual stocks appreciate explosively relative to the market.

### 3. Fama-French Three-Factor Idiosyncratic Component

Removes systematic exposure to:

- Market risk
- Size (SMB)
- Value (HML)

The accumulated regression residual represents the idiosyncratic component of the stock.

### 4. Fama-French Five-Factor Idiosyncratic Component

Extends the previous specification by additionally controlling for:

- Profitability (RMW)
- Investment (CMA)

This provides the strictest test of firm-specific speculative exuberance used in the project.

---

## Econometric Methodology

### Bubble Detection

The project applies the recursive right-tailed unit-root methodology developed by Phillips, Shi and Yu (2015).

The PSY/GSADF framework is designed to identify periodically collapsing explosive episodes that conventional full-sample ADF tests may fail to detect.

### Bubble Severity

Following Phillips (2023), the analysis also examines the localizing-rate parameter associated with near-unit-root processes.

The estimated dynamics are interpreted in terms of:

- Mildly Integrated
- Local-to-Unity
- Mildly Explosive

This allows the analysis to go beyond binary bubble detection and examine the magnitude of explosive behavior.

---

## Sample Design

The analysis is divided into two periods.

Period 1: January 2010 – November 2021

The PSY real-time detection methodology is used to identify and date speculative episodes.

Period 2: December 2022 – June 2026

The GSADF statistic is used to test for explosive dynamics over the shorter post-2022 sample.

---

## Key Findings

The results provide the strongest evidence of speculative exuberance for Microsoft and NVIDIA during Period 1.

Bubble signals for Apple, Tesla and Meta are more dependent on the specification used, while Alphabet and Amazon generally show weaker evidence after controlling for systematic factors.

For Period 2, none of the Magnificent Seven produces statistically significant GSADF evidence of explosive behavior across the four specifications.

Near-unit-root severity analysis generally indicates strong persistence rather than statistically confirmed Mildly Explosive behavior.

---

## Main References

Phillips, P. C. B., Shi, S., & Yu, J. (2015). Testing for multiple bubbles: Historical episodes of exuberance and collapse in the S&P 500.

Phillips, P. C. B. (2023). Estimation and inference with near unit roots.

Basele, R. B., Phillips, P. C. B., & Shi, S. (2025). Speculative Bubbles in the Recent AI Boom: Nasdaq and the Magnificent Seven.

---

## Software

The empirical analysis was conducted in R.

Main packages include:

- psymonitor
- exuber
- tidyquant
- frenchdata
- ggplot2
- dplyr

---

## Repository Structure

```text
R/
└── stocks/
    ├── 01_AAPL_analysis.R
    ├── 02_MSFT_analysis.R
    ├── 03_GOOGL_analysis.R
    ├── 04_AMZN_analysis.R
    ├── 05_NVDA_analysis.R
    ├── 06_META_analysis.R
    └── 07_TSLA_analysis.R

data/
├── revenue/
└── market/

results/
├── figures/
└── tables/

thesis/
└── Asset_Bubble_Detection_Thesis.pdf
```

---

## Author

Armin Chiani  
Master's Degree in Economics  
University of Turin
