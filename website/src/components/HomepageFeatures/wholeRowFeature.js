import React from 'react';
import styles from './styles.module.css';
import useBaseUrl from '@docusaurus/useBaseUrl';
import clsx from 'clsx';
import Link from '@docusaurus/Link';

function WholeRowFeature({ imageUrl, title, description, buttons, contentAlignment, imageAlignment, colSize = 'col--6', }) {
  const resolvedImageUrl = useBaseUrl(imageUrl);

  const imageAlignmentClass =
    imageAlignment === 'right' ? styles.alignRight :
    imageAlignment === 'center' ? styles.alignCenter :
    styles.alignLeft;

  const contentAlignmentClass =
    contentAlignment === 'right' ? styles.alignRight :
    contentAlignment === 'center' ? styles.alignCenter :
    styles.alignLeft;

  return (
    <div className={clsx('col', colSize, styles.featureContainer)}>
      {resolvedImageUrl && (
        <div className={clsx(styles.featureImage, imageAlignmentClass)}>
          <img
            className={styles.largeFeatureImage}
            src={resolvedImageUrl}
            alt={title || ''}
          />
        </div>
      )}

      <div className={clsx(styles.featureContent, contentAlignmentClass)}>
        {title && <h3>{title}</h3>}
        {description && <p>{description}</p>}
        {buttons && buttons.length > 0 && (
          <div className={clsx(styles.buttonContainer, contentAlignmentClass)}>
            {buttons.map((button, index) => (
              <Link key={index} className="button button--primary" to={button.buttonLink}>
                {button.buttonText}
              </Link>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}

export default WholeRowFeature;
