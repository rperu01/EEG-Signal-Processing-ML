# EEG-Signal-Processing-ML
EEG-based ML project focused on detecting mental arithmetic tasks from noisy physiological signals. Preprocessed EEG data using band-pass filtering and FFT, extracted frequency-domain features (beta band power, statistics), and trained an SVM classifier, achieving 61% accuracy on a public dataset.

## Project Overview
This project aimed to detect mental arithmetic tasks from EEG signals using machine learning. I worked with a publicly available EEG dataset, preprocessed low-SNR signals, and developed a classifier to identify when participants were performing mental arithmetic based on their brain activity patterns.

## My Role
I led the signal processing and machine learning pipeline:

Preprocessed EEG signals using MATLAB, applied band-pass filters, and performed FFT to extract frequency-domain features.

Extracted features such as beta-band power, standard deviation, and task labels.

Built and trained a Support Vector Machine (SVM) classifier, and analyzed model performance and accuracy.

## Methods & Tools
Signal Processing: MATLAB, FFT, band-pass filtering (13–30 Hz beta waves)

Feature Extraction: Relative band power, standard deviation, task labels, participant metadata

Machine Learning: SVM classifier, cross-validation, accuracy assessment

Data: Public EEG dataset, 61% classification accuracy achievedSignal Processing: MATLAB, FFT, band-pass filtering (13–30 Hz beta waves)

Feature Extraction: Relative band power, standard deviation, task labels, participant metadata

Machine Learning: SVM classifier, cross-validation, accuracy assessment

Data: Public EEG dataset, 61% classification accuracy achieved

## Results
Successfully detected mental arithmetic tasks from EEG signals with 61% accuracy.

Demonstrated ability to process low-SNR physiological data and extract meaningful features.

Pipeline showed potential for extending to other cognitive-state detection applications.

## Notes
Data was publicly available; no confidential information was used.

Focused on reproducible methods and clear documentation to showcase signal processing and ML workflow.

Igor Zyma, Ivan Seleznov, Anton Popov, Mariia Chernykh, and Oleksii
Shpenkov, “Eeg during mental arithmetic tasks,” https://physionet.org/
content/eegmat/1.0.0/, 2018, accessed: Oct, 12 2024.
