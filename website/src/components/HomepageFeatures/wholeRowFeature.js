import React from 'react';
import styles from './styles.module.css';
import useBaseUrl from '@docusaurus/useBaseUrl';
import clsx from 'clsx';
import Link from '@docusaurus/Link';

function WholeRowFeature({ imageUrl, title, description, buttons, contentAlignment, imageAlignment }) {
  const resolvedImageUrl = useBaseUrl(imageUrl);

  return (
    <div className={clsx('col col--6', styles.featureContainer, imageAlignment === 'right' ? styles.alignRight : styles.alignLeft)}>
      {resolvedImageUrl && (
        <div className={styles.featureImage}>
          <img
            className={styles.largeFeatureImage}
            src={resolvedImageUrl}
            alt={title || ''}
          />
        </div>
      )}

      <div className={clsx(styles.featureContent, contentAlignment === 'right' ? styles.alignRight : '')}>
        {title && <h3>{title}</h3>}
        {description && <p>{description}</p>}
        {buttons && buttons.length > 0 && (
          <div className={clsx(styles.buttonContainer, contentAlignment === 'right' ? styles.alignRight : '')}>
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
