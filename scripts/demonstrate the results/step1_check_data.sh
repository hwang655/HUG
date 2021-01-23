#!/bin/sh

R CMD BATCH --no-save --no-restore "--args ind='5278'" step1_check_data.R step1_check_data_5278_ASD.Rout &

R CMD BATCH --no-save --no-restore "--args ind='5531'" step1_check_data.R step1_check_data_5531_ASD.Rout &

R CMD BATCH --no-save --no-restore "--args ind='4341'" step1_check_data.R step1_check_data_4341_control.Rout &

R CMD BATCH --no-save --no-restore "--args ind='5958'" step1_check_data.R step1_check_data_5958_control.Rout &
