import React from 'react';
import styles from './styles.module.css';
import useBaseUrl from '@docusaurus/useBaseUrl';

function FeatureGallery({ header, images, description }) {
  return (
    <div className={styles.featureGallery}>
      <h2 className={styles.header}>{header}</h2>

      <div className={styles.imageGallery}>
        {images.map((image, index) => {
          const resolvedImageUrl = useBaseUrl(image.src);

          return (
            <a key={index} href={image.link} target="_blank" rel="noopener noreferrer">
              <img className={styles.galleryImage} src={resolvedImageUrl} alt={image.alt || `Image ${index + 1}`} />
            </a>
          );
        })}
      </div>

      {description && <p className={styles.description}>{description}</p>}
    </div>
  );
}

export default FeatureGallery;
